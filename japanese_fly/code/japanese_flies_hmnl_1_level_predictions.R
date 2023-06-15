library(rstan)
library(MASS)
library(tidyverse)
library(ggtern)
library(here)

source(here("japanese_fly/code/japanese_flies_hmnl_1_level_utils.R"))

# Update to correct one
# model_stan = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_4_stan_object_30images.rds"))
model_stan = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_7_stan_object_30images.rds"))

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




betas_level_0_summary %>% 
  mutate(exp = "All experiments", line_thickness_1 = 2, line_thickness_2 = 1) %>% 
  bind_rows(betas_level_1_summary %>% 
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





betas_level_1_summary %>% 
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



betas_level_1_summary %>% 
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



predict_utilities_design = function(
  df_to_predict, 
  beta_level_1_draws_hmnl = NULL,
  beta_level_0_draws_hmnl = NULL, 
  Sigma_level_0_draws_hmnl = NULL,
  n_draws_per_posterior_sample = NULL,
  include_all_utilities = T, 
  seed = NULL
){
  
  # Example of how to use by generating the beta_level_1 on the fly
  # predicted_utilities_df_01 = predict_utilities_design(
  #   df_to_predict = df_to_predict_01,
  #   beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  #   Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  #   include_all_utilities = F,
  #   n_draws_per_posterior_sample = 10,
  #   seed = 2023
  # )
  
  # Example of how to use with a matrix of precomputed betas_level_1
  # beta_level_1_draws_hmnl_10samples = create_level_1_betas_hmnl(
  #   beta_level_0_draws_hmnl, Sigma_level_0_draws_hmnl, n_draws_per_posterior_sample = 10,
  #   seed = 2023
  # )
  # 
  # predicted_utilities_df_01 = predict_utilities_design(
  #   df_to_predict = df_to_predict_01, 
  #   beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  #   include_all_utilities = F
  # )
  
  X_to_predict_list <- create_model_matrix_second_order_scheffe(df_to_predict)
  
  if(!is.null(beta_level_1_draws_hmnl)){
    
    if(!is.null(beta_level_0_draws_hmnl)) warning("beta_level_0_draws_hmnl provided but ignored")
    if(!is.null(Sigma_level_0_draws_hmnl)) warning("Sigma_level_0_draws_hmnl provided but ignored")
    if(!is.null(n_draws_per_posterior_sample)) warning("n_draws_per_posterior_sample provided but ignored")
    
    n_draws = nrow(beta_level_1_draws_hmnl)
    utilities_post = matrix(0.0, ncol = n_draws, nrow = nrow(X_to_predict_list$X))
    
    for(d in 1:n_draws){
      beta_level_1 = beta_level_1_draws_hmnl[d,]
      utilities_d = X_to_predict_list$X %*% beta_level_1
      utilities_post[, d] = utilities_d
    }
    
    
    
  } else{
    
    n_draws = nrow(beta_level_0_draws_hmnl[[1]])
    utilities_post = matrix(0.0, ncol = n_draws_per_posterior_sample*n_draws, nrow = nrow(X_to_predict_list$X))
    
    counter = 1
    if(!is.null(seed)) set.seed(seed = seed)
    for(d in 1:n_draws){
      Sigma_level_0_d = Sigma_level_0_draws_hmnl[[1]][d, , ]
      beta_level_0_d = beta_level_0_draws_hmnl[[1]][d, ]
      beta_level_1 = mvrnorm(n_draws_per_posterior_sample, beta_level_0_d, Sigma_level_0_d)
      
      if(n_draws_per_posterior_sample > 1){
        for(i in 1:n_draws_per_posterior_sample){
          utilities_d = X_to_predict_list$X %*% beta_level_1[i, ]
        }
      } else{
        utilities_d = X_to_predict_list$X %*% beta_level_1
      }
      
      
      
      utilities_post[, counter:(counter + n_draws_per_posterior_sample - 1)] = utilities_d
      counter = counter + n_draws_per_posterior_sample
    }
    
  } # end if
  
  
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
  
  return(out)
  
}








beta_level_1_draws_hmnl_10samples = create_level_1_betas_hmnl(
  beta_level_0_draws_hmnl, Sigma_level_0_draws_hmnl, n_draws_per_posterior_sample = 10,
  seed = 2023
)








lattice_design_7_ing_01 = create_lattice_design(n_var = 7, n_levels = 7) %>% 
  set_names(c('R','O','Y','G','B','P', 'UV'))


df_to_predict_01 = expand_grid(
  lattice_design_7_ing_01, 
  intensity = seq(-1, 1, by = 0.2),
  no_choice = 0
) 


predicted_utilities_df_01 = predict_utilities_design(
  df_to_predict = df_to_predict_01, 
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = F
)



predicted_utilities_df_01$utilities_summary %>% 
  arrange(desc(p50))



predicted_utilities_df_01$utilities_summary %>% 
  top_n(100, p50) %>% 
  arrange(desc(p50)) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()


# Sample 10000 points to see utilities and posterior intervals
predicted_utilities_df_01$utilities_summary %>% 
  sample_n(10000) %>% 
  arrange(desc(p50)) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point(size = 0.05) +
  geom_linerange(aes(ymin = p10, ymax = p90), linewidth = 0.1, alpha = 0.6) +
  theme_bw()


# Best vs worst
predicted_utilities_df_01$utilities_summary %>% 
  top_n(500, p50) %>% 
  bind_rows(
    predicted_utilities_df_01$utilities_summary %>% 
      top_n(-500, p50)
  ) %>% 
  arrange(desc(p50)) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point(size = 0.1) +
  geom_linerange(aes(ymin = p10, ymax = p90), linewidth = 0.4, alpha = 0.7) +
  theme_bw()




# Top 20 mixture values
top_20_mixtures_01 = predicted_utilities_df_01$utilities_summary %>% 
  arrange(desc(p50)) %>% 
  select(R, O, Y, G, B, P, UV) %>% 
  distinct() %>% 
  head(20) %>% 
  as.data.frame() 


top_20_mixtures_01_min = top_20_mixtures_01 %>% 
  apply(., 2, min)


top_20_mixtures_01_max = top_20_mixtures_01 %>% 
  apply(., 2, max)


top_20_mixtures_01

# Range of values of the mixtures with highest utilities
top_20_mixtures_01_min
top_20_mixtures_01_max


# More granular around the best values predicted in previous step
lattice_design_7_ing_02 = create_lattice_design(
  n_var = 7, n_levels = 7, 
  limits = rbind(top_20_mixtures_01_min, top_20_mixtures_01_max)
) %>% 
  set_names(c('R','O','Y','G','B','P', 'UV'))






df_to_predict_02 = expand_grid(
  lattice_design_7_ing_02, 
  intensity = seq(-1, 1, by = 0.1),
  no_choice = 0
) 


# 78 seconds with 600 samples, 10 draws per sample, and 67K rows in df_to_predict
# 20 seconds with 1200 samples, 10 draws per sample, and 11K rows in df_to_predict
Sys.time()
predicted_utilities_df_02 = predict_utilities_design(
  df_to_predict = df_to_predict_02, 
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = F
)
Sys.time()


predicted_utilities_df_02$utilities_summary %>% 
  arrange(desc(p50))


predicted_utilities_df_02$utilities_summary %>% 
  top_n(100, p50) %>% 
  arrange(desc(p50)) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()


# Sample 6000 points to see utilities and posterior intervals
predicted_utilities_df_02$utilities_summary %>% 
  sample_n(6000) %>% 
  arrange(desc(p50)) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point(size = 0.05) +
  geom_linerange(aes(ymin = p10, ymax = p90), linewidth = 0.2, alpha = 0.7) +
  theme_bw()


# Best vs worst
predicted_utilities_df_02$utilities_summary %>% 
  top_n(500, p50) %>% 
  bind_rows(
    predicted_utilities_df_02$utilities_summary %>% 
      top_n(-500, p50)
  ) %>% 
  arrange(desc(p50)) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point(size = 0.1) +
  geom_linerange(aes(ymin = p10, ymax = p90), linewidth = 0.4, alpha = 0.7) +
  theme_bw()



# Top 20 mixture values
top_20_mixtures_02 = predicted_utilities_df_02$utilities_summary %>% 
  arrange(desc(p50)) %>% 
  select(R, O, Y, G, B, P, UV) %>% 
  distinct() %>% 
  head(20) %>% 
  as.data.frame() 


top_20_mixtures_02_min = top_20_mixtures_02 %>% 
  apply(., 2, min)


top_20_mixtures_02_max = top_20_mixtures_02 %>% 
  apply(., 2, max)


top_20_mixtures_02

# Range of values of the mixtures with highest utilities
top_20_mixtures_02_min
top_20_mixtures_02_max


# More granular around the best values predicted in previous step
lattice_design_7_ing_02 = create_lattice_design(
  n_var = 7, n_levels = 8, 
  limits = rbind(top_20_mixtures_02_min, top_20_mixtures_02_max)
) %>% 
  set_names(c('R','O','Y','G','B','P', 'UV'))



df_to_predict_03 = lattice_design_7_ing_02 %>% 
  expand_grid(
    .,
    intensity = seq(-1, 1, by = 0.05), 
    no_choice = 0
  )




# 34 seconds with 600 samples, 10 draws per sample, and 20K rows in df_to_predict
# 90 seconds with 1200 samples, 10 draws per sample, and 37K rows in df_to_predict
Sys.time()
predicted_utilities_df_03 = predict_utilities_design(
  df_to_predict = df_to_predict_03, 
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = F
)
Sys.time()


predicted_utilities_df_03$utilities_summary %>% 
  arrange(desc(p50))





predicted_utilities_df_03$utilities_summary %>% 
  top_n(100, p50) %>% 
  arrange(desc(p50)) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()



predicted_utilities_df_03$utilities_summary %>% 
  top_n(200, p50) %>% 
  bind_rows(
    predicted_utilities_df_03$utilities_summary %>% 
      top_n(-200, p50)
  ) %>% 
  arrange(desc(p50)) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()



# Top 20 mixture values
top_20_mixtures_03 = predicted_utilities_df_03$utilities_summary %>% 
  arrange(desc(p50)) %>% 
  select(R, O, Y, G, B, P, UV) %>% 
  distinct() %>% 
  head(20) %>% 
  as.data.frame() 


top_20_mixtures_03_min = top_20_mixtures_03 %>% 
  apply(., 2, min)


top_20_mixtures_03_max = top_20_mixtures_03 %>% 
  apply(., 2, max)


top_20_mixtures_03

# Range of values of the mixtures with highest utilities
top_20_mixtures_03_min
top_20_mixtures_03_max







# top 5 overall
top_5_03_overall = predicted_utilities_df_03$utilities_summary %>% 
  top_n(5, p50)


# top 5 with lower SD
top_5_03_low_sd = predicted_utilities_df_03$utilities_summary %>% 
  top_n(100, p50) %>% 
  arrange(sd) %>% 
  head(5)



top_5_03_overall %>% 
  bind_rows(
    top_5_03_low_sd
  ) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()







top_5_03_min = top_5_03_low_sd %>% 
  bind_rows(top_5_03_overall) %>% 
  select(R, O, Y, G, B, P, UV, intensity) %>% 
  apply(., 2, min)


top_5_03_max = top_5_03_low_sd %>% 
  bind_rows(top_5_03_overall) %>% 
  select(R, O, Y, G, B, P, UV, intensity) %>% 
  apply(., 2, max)


top_20_mixtures_03


top_5_03_min
top_5_03_max








# Even more granular around the best values predicted in previous step
lattice_design_7_ing_03 = create_lattice_design(
  n_var = 7, n_levels = 10, 
  limits = as.matrix(bind_rows(top_5_03_min, top_5_03_max) %>% select(-intensity))
) %>% 
  set_names(c('R','O','Y','G','B','P', 'UV'))



df_to_predict_04 = lattice_design_7_ing_03 %>% 
  expand_grid(
    .,
    intensity = seq(top_5_03_min["intensity"], top_5_03_max["intensity"], by = 0.01), 
    no_choice = 0
  )




# 32 seconds with 600 samples, 10 draws per sample, and 28K rows in df_to_predict
# 25 seconds with 1200 samples, 10 draws per sample, and 11K rows in df_to_predict
Sys.time()
predicted_utilities_df_04 = predict_utilities_design(
  df_to_predict = df_to_predict_04, 
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = F
)
Sys.time()


predicted_utilities_df_04$utilities_summary %>% 
  arrange(desc(p50))








x_vector = c(0.07, 0.2, 0, 0, 0, 0.23, 0.5, 0.65)


get_utility_distribution = function(
  x_vector, beta_level_1_draws_hmnl, include_all_utilities = T
){
  
  df_to_predict = as.data.frame(t(x_vector)) %>% 
    set_names(c("R", "O", "Y", "G", "B", "P", "UV", "intensity")) %>% 
    mutate(no_choice = 0)
  
  pred_list = predict_utilities_design(
    df_to_predict, 
    beta_level_1_draws_hmnl = beta_level_1_draws_hmnl,
    include_all_utilities = include_all_utilities)
  
  return(pred_list)
  
}



get_utility_distribution(
  x_vector = c(0.1, 0.114, 0, 0, 0, 0.3, 0.486, 0.65),
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = F
)


get_utility_distribution(
  x_vector = c(0.07, 0.2, 0, 0, 0, 0.23, 0.5, 0.65),
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = F
)














utility_optim_function = function(x_vector_optim){
  # First element is 1 minus the sum of the rest
  n_param = length(x_vector_optim)
  
  x_vector2 = rep(0.0, n_param+1)
  x_vector2[1] = 1 - sum(x_vector_optim[1:(n_param-1)])
  x_vector2[2:(n_param+1)] = x_vector_optim
  
  return(
    get_utility_distribution(
      x_vector = x_vector2,
      beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
      include_all_utilities = F
    )$utilities_summary$p50
  )
}

utility_optim_function(c(0.114, 0, 0, 0, 0.3, 0.486, 0.65))

transform_optim_to_mixture_pv = function(param){
  n_param = length(param)
  opt_mixt = rep(0.0, n_param+1)
  opt_mixt[1] = 1 - sum(param[1:(n_param-1)])
  opt_mixt[2:(n_param+1)] = param
  return(opt_mixt)
}




Sys.time()
# 2 minutes
optim_mixt_result_01 = optim(
  par = c(0.2, 0, 0, 0, 0.2, 0.5, 0.65),
  fn = utility_optim_function,
  lower = c(0, 0, 0, 0, 0, 0, -1),
  upper = c(1, 1, 1, 1, 1, 1, 1),
  method = "L-BFGS-B",
  control = list(fnscale = -1)
)
Sys.time()
optim_mixt_result_01


round(transform_optim_to_mixture_pv(optim_mixt_result_01$par), 4)

get_utility_distribution(
  x_vector = transform_optim_to_mixture_pv(optim_mixt_result_01$par),
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = F
)






optim_initial_par_02 = predicted_utilities_df_04$utilities_summary %>% 
  top_n(1, p50) %>% 
  select(O, Y, G, B, P, UV, intensity) %>% 
  as.numeric()



Sys.time()
# 15 seconds
optim_mixt_result_02 = optim(
  par = optim_initial_par_02,
  fn = utility_optim_function,
  lower = c(0, 0, 0, 0, 0, 0, -1),
  upper = c(1, 1, 1, 1, 1, 1, 1),
  method = "L-BFGS-B",
  control = list(fnscale = -1)
)
Sys.time()
optim_mixt_result_02


round(transform_optim_to_mixture_pv(optim_initial_par_02), 4)
round(transform_optim_to_mixture_pv(optim_mixt_result_02$par), 4)

get_utility_distribution(
  x_vector = transform_optim_to_mixture_pv(optim_mixt_result_02$par),
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = F
)









optim_initial_par_03 = rep(0, 7)

Sys.time()
# 1 minute
optim_mixt_result_03 = optim(
  par = optim_initial_par_03,
  fn = utility_optim_function,
  lower = c(0, 0, 0, 0, 0, 0, -1),
  upper = c(1, 1, 1, 1, 1, 1, 1),
  method = "L-BFGS-B",
  control = list(fnscale = -1)
)
Sys.time()
optim_mixt_result_03


round(transform_optim_to_mixture_pv(optim_initial_par_03), 4)
round(transform_optim_to_mixture_pv(optim_mixt_result_03$par), 4)

get_utility_distribution(
  x_vector = transform_optim_to_mixture_pv(optim_mixt_result_03$par),
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = F
)




best_utilities_so_far = predicted_utilities_df_01$utilities_summary %>% 
  top_n(5, p50) %>% 
  bind_rows(
    predicted_utilities_df_02$utilities_summary %>% 
      top_n(5, p50)
    
  ) %>% 
  bind_rows(
    predicted_utilities_df_03$utilities_summary %>% 
      top_n(5, p50)
    
  ) %>% 
  bind_rows(
    predicted_utilities_df_04$utilities_summary %>% 
      top_n(5, p50)
    
  ) %>% 
  bind_rows(
    get_utility_distribution(
      x_vector = transform_optim_to_mixture_pv(optim_mixt_result_01$par),
      beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
      include_all_utilities = F
    )$utilities_summary
  ) %>% 
  bind_rows(
    get_utility_distribution(
      x_vector = transform_optim_to_mixture_pv(optim_mixt_result_02$par),
      beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
      include_all_utilities = F
    )$utilities_summary
  ) %>% 
  bind_rows(
    get_utility_distribution(
      x_vector = transform_optim_to_mixture_pv(optim_mixt_result_03$par),
      beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
      include_all_utilities = F
    )$utilities_summary
  ) %>% 
  arrange(desc(p50))

best_utilities_so_far

best_utilities_so_far %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()





# # Probably not necessary anymore
# get_utility_plugin_distribution = function(
#   x_vector, beta_level_0_draws_hmnl, Sigma_level_0_draws_hmnl, include_all_utilities = T
# ){
#   
#   df_to_predict = as.data.frame(t(x_vector)) %>% 
#     set_names(c("R", "O", "Y", "G", "B", "P", "UV", "intensity")) %>% 
#     mutate(no_choice = 0)
#   
#   X_to_predict_list <- create_model_matrix_second_order_scheffe(df_to_predict)
#   
#   n_draws = nrow(beta_level_0_draws_hmnl[[1]])
#   
#   utilities_post = matrix(0.0, ncol = n_draws, nrow = nrow(X_to_predict_list$X))
#   
#   for(d in 1:n_draws){
#     beta_level_0_d = beta_level_0_draws_hmnl[[1]][d, ]
#     utilities_d = X_to_predict_list$X %*% beta_level_0_d
#     
#     utilities_post[, d] = utilities_d
#   }
#   
#   
#   quantiles_utilities = t(apply(utilities_post, 1, function(x) c(quantile(x, probs = c(0.1, 0.5, 0.9)), sd = sd(x)) ))
#   
#   
#   out_utilities_summary = df_to_predict %>% 
#     mutate(order_original = 1:n()) %>% 
#     bind_cols(quantiles_utilities %>% 
#                 as_tibble() %>% 
#                 set_names(c("p10", "p50", "p90", "sd")))
#   
#   
#   if(include_all_utilities){
#     out = list(
#       utilities_summary = out_utilities_summary,
#       utilities_post = utilities_post
#     )
#   }else{
#     out = list(
#       utilities_summary = out_utilities_summary
#     )
#   }
#   
#   return(out)
#   
# }
# 
# 
# 
# 
# 
# get_utility_plugin_distribution(
#   x_vector = c(0.1, 0.114, 0, 0, 0, 0.3, 0.486, 0.65),
#   beta_level_0_draws_hmnl = beta_level_0_draws_hmnl, 
#   Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
#   include_all_utilities = F
# )
# 
# 
# 
# plugin_utility_optim_function = function(x_vector_optim){
#   # First element is 1 minus the sum of the rest
#   n_param = length(x_vector_optim)
#   
#   x_vector2 = rep(0.0, n_param+1)
#   x_vector2[1] = 1 - sum(x_vector_optim[1:(n_param-1)])
#   x_vector2[2:(n_param+1)] = x_vector_optim
#   
#   return(
#     get_utility_plugin_distribution(
#       x_vector = x_vector2,
#       beta_level_0_draws_hmnl = beta_level_0_draws_hmnl, 
#       Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
#       include_all_utilities = F
#     )$utilities_summary$p50
#   )
# }
# 
# plugin_utility_optim_function(c(0.1, 0.114, 0, 0, 0, 0.3, 0.486, 0.65))
