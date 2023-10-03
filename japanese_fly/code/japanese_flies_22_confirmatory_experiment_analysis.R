library(rstan)
library(MASS)
library(tidyverse)
library(ggtern)
library(here)

source(here("japanese_fly/code/japanese_flies_hmnl_1_level_utils.R"))

counts_japanese = read_csv(here("japanese_fly/out/counts_japanese_choice_format.csv"))

model_stan = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_9_stan_object_30images.rds"))

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
  mutate(variable = fct_reorder(variable, ix, .desc = F))


########################################
##### Data plots
########################################

hand_picked_colors = readRDS(here("japanese_fly/out/handpicked_colors_5th_i_opt_design.rds"))

design_confirmatory_ordered = read.csv(here("japanese_fly/out/confirmatory_design_japanese_ordered.csv"))

design_confirmatory_with_mixture = read.csv(here("japanese_fly/out/confirmatory_design_japanese_mapping.csv"))

# The first 60 rows of design_confirmatory_ordered have the choice sets of interest,
# but the rows were randomized, so cs_new has the correct rows of the randomized design
choice_sets_of_interest = unique(design_confirmatory_ordered$cs_new[1:60])


# # The first 22 rows here are the 22 choice sets that I created and forcefully added in the experiment
# # The last 6 rows have one fixed alternative of the 6 handpicked colors I wanted to check out
# hand_picked_choice_sets = design_confirmatory_with_mixture[choice_sets_of_interest,]






choice_sets_to_analyze = design_confirmatory_with_mixture[choice_sets_of_interest, "folder"]
# choice_sets_to_analyze = design_confirmatory_with_mixture[1:28, "folder"]





filtered_exp_10_cs_of_interest = counts_japanese %>% 
  filter(experiment == 10, folder %in% choice_sets_to_analyze) %>% 
  select(R:no_choice)

filtered_exp_10_cs_of_interest_to_predict_list <- create_model_matrix_second_order_scheffe(filtered_exp_10_cs_of_interest)


predicted_point_utilities = filtered_exp_10_cs_of_interest_to_predict_list$X %*% betas_level_0_summary$mean


filtered_exp_10_cs_of_interest_utilities = counts_japanese %>% 
  filter(experiment == 10, folder %in% choice_sets_to_analyze) %>% 
  mutate(pred_ut = as.numeric(predicted_point_utilities))



filtered_exp_10_cs_of_interest_utilities

filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") 




choice_sets_of_interest_catalog = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  select(folder, choice_set, alternative, pred_ut, R:intensity) %>% 
  distinct() %>% 
  mutate(color = paste(R, O, Y, G, B, P, UV, intensity)) %>% 
  mutate(color_ix = as.integer(as.factor(color)))


exp_10_handpicked_colors_indices_and_utilities = choice_sets_of_interest_catalog %>% 
  mutate(exp_pred_ut = exp(pred_ut)) %>% 
  select(color_ix, pred_ut, exp_pred_ut, R:intensity) %>% 
  distinct() %>% 
  arrange(color_ix)

exp_10_handpicked_colors_indices_and_utilities

write_csv(exp_10_handpicked_colors_indices_and_utilities,
          here("japanese_fly/out/plots/exp_10_handpicked_colors_indices_and_utilities.csv"))



width_counts_plots = 35
height_counts_plots = 23


counts_exp_10_handpicked_colors_01 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  ggplot() +
  geom_jitter(aes(x = alternative, y = count), size = 0.6) +
  facet_wrap(~folder) +
  theme_bw()


ggsave(
  plot = counts_exp_10_handpicked_colors_01,
  filename = here("japanese_fly/out/plots/counts_exp_10_handpicked_colors_01.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)




utilities_counts_tibble = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  group_by(folder, choice_set, alternative) %>% 
  summarize(
    min_counts = min(count),
    max_counts = max(count),
    p10_counts = quantile(count, 0.1),
    p90_counts = quantile(count, 0.9),
    median_counts = median(count),
    pred_ut = unique(pred_ut)) %>%
  # if ratio_ut_l_by_r > 1 then left is more preferred
  mutate(ratio_ut_l_by_r = exp(first(pred_ut) - last(pred_ut)),
         LU = round(exp(first(pred_ut)), 2),
         RU = round(exp(last(pred_ut)), 2)) %>% 
  ungroup() %>% 
  arrange(desc(ratio_ut_l_by_r), folder) %>% 
  mutate(cs_ordered = rep(1:(nrow(.)/2), each = 2)) %>% 
  # mutate(group = paste0(sprintf("%02d", cs_ordered), ") ", round(ratio_ut_r_by_l, 2))) %>%
  mutate(group = paste0(sprintf("%02d", cs_ordered), ") LU=", LU,
                        " / RU=", RU, " (ratio=", round(ratio_ut_l_by_r, 2), ")"))



counts_exp_10_handpicked_colors_02 = utilities_counts_tibble %>% 
  ggplot(aes(x = alternative, y = median_counts)) +
  geom_point() +
  geom_linerange(aes(ymin = p10_counts, ymax = p90_counts)) +
  facet_wrap(~group) +
  theme_bw() +
  ggtitle("Percentiles 10 and 90 of counts in each of the choice sets of interest")


ggsave(
  plot = counts_exp_10_handpicked_colors_02,
  filename = here("japanese_fly/out/plots/counts_exp_10_handpicked_colors_02.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)



counts_exp_10_handpicked_colors_03 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  left_join(
    utilities_counts_tibble %>% 
      select(folder, alternative, group)
  ) %>% 
  ggplot() +
  geom_jitter(aes(x = alternative, y = count), size = 0.3) +
  facet_wrap(~group) +
  theme_bw()

ggsave(
  plot = counts_exp_10_handpicked_colors_03,
  filename = here("japanese_fly/out/plots/counts_exp_10_handpicked_colors_03.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)



counts_exp_10_handpicked_colors_04 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  left_join(
    utilities_counts_tibble %>% 
      select(folder, alternative, group)
  ) %>% 
  left_join(
    choice_sets_of_interest_catalog %>% 
      select(folder, alternative, color_ix)
  ) %>% 
  mutate(x = paste0(substr(alternative, 1, 1), ") color: ", color_ix)) %>% 
  ggplot() +
  geom_jitter(aes(x = x, y = count), size = 0.3) +
  facet_wrap(~group, scales = "free_x") +
  xlab("Alternative") +
  theme_bw() +
  xlab("") +
  ylab("")

ggsave(
  plot = counts_exp_10_handpicked_colors_04,
  filename = here("japanese_fly/out/plots/counts_exp_10_handpicked_colors_04.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)




counts_exp_10_handpicked_colors_05 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  left_join(
    utilities_counts_tibble %>% 
      select(folder, alternative, group, min_counts, max_counts, median_counts)
  ) %>% 
  left_join(
    choice_sets_of_interest_catalog %>% 
      select(folder, alternative, color_ix)
  ) %>% 
  mutate(x = paste0(substr(alternative, 1, 1), ") color: ", color_ix)) %>% 
  ggplot(aes(x = x)) +
  geom_jitter(aes(y = count), size = 0.2, width = 0.3, alpha = 0.8) +
  geom_point(aes(y = median_counts), color = "red", size = 2.3) +
  geom_linerange(aes(ymin = min_counts, ymax = max_counts), color = "red", size = 1.3) +
  facet_wrap(~group, scales = "free_x") +
  xlab("Alternative") +
  theme_bw() +
  xlab("") +
  ylab("")


ggsave(
  plot = counts_exp_10_handpicked_colors_05,
  filename = here("japanese_fly/out/plots/counts_exp_10_handpicked_colors_05.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)





counts_exp_10_handpicked_colors_06 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  left_join(
    utilities_counts_tibble %>% 
      select(folder, alternative, group, min_counts, max_counts, median_counts)
  ) %>% 
  left_join(
    choice_sets_of_interest_catalog %>% 
      select(folder, alternative, color_ix)
  ) %>% 
  mutate(x = paste0(substr(alternative, 1, 1), ") color: ", color_ix)) %>% 
  ggplot(aes(x = x)) +
  geom_jitter(aes(y = count), size = 0.2, width = 0.3, alpha = 0.8) +
  geom_point(aes(y = median_counts), color = "red", size = 2.3) +
  geom_linerange(aes(ymin = min_counts, ymax = max_counts), color = "red", size = 1.3) +
  geom_text(aes(y = 1.1*max(max_counts), label = paste0("Median = ", median_counts)),
            check_overlap = T, vjust = 1.1) +
  facet_wrap(~group, scales = "free_x") +
  theme_bw() +
  xlab("") +
  ylab("")


ggsave(
  plot = counts_exp_10_handpicked_colors_06,
  filename = here("japanese_fly/out/plots/counts_exp_10_handpicked_colors_06.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)








### Similar plot with predicted counts

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
        utilities_d = X_to_predict_list$X %*% t(beta_level_1)
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



predict_N_flies_design = function(
  choice_sets,
  utilities_post,
  n_flies = 80,
  percentiles = c(10, 25, 50, 75, 90),
  include_all_simulated_responses = T,
  seed = NULL
){
  
  J = length(choice_sets)/length(unique(choice_sets))
  
  n_rows = nrow(utilities_post)
  n_draws = ncol(utilities_post)
  
  
  probs_post = array(0.0, c(n_draws, J, n_rows/J))
  N_flies_post = array(0.0, dim(probs_post))
  
  if(is.null(seed)) set.seed(seed)
  
  for(d in 1:n_draws){
    
    probs_d = sapply(split(utilities_post[, d], choice_sets), function(x){
      exp_x = exp(x)
      exp_x/sum(exp_x)
    })
    
    for(i in 1:ncol(probs_d)){
      N_flies_post[d, , i] = as.numeric(rmultinom(n = 1, size = n_flies, prob = probs_d[, i]))
    }
  }
  
  
  mean_N_flies_post = apply(N_flies_post, c(2, 3), mean)
  # median_N_flies_post = apply(N_flies_post, c(2, 3), median)
  # p10_N_flies_post = apply(N_flies_post, c(2, 3), quantile, 0.1)
  # p90_N_flies_post = apply(N_flies_post, c(2, 3), quantile, 0.9)
  
  N_flies_post_percentiles = apply(N_flies_post, c(2, 3), quantile, percentiles/100)
  
  N_flies_post_percentiles_df = matrix(rep(0.0, length(N_flies_post_percentiles)), ncol = length(percentiles)) %>% 
    as.data.frame() %>% 
    set_names(paste0("p", percentiles, "_N_flies"))
  
  for(i in 1:ncol(N_flies_post_percentiles_df)){
    N_flies_post_percentiles_df[, i] = as.numeric(N_flies_post_percentiles[i, , ])
  }
  
  
  N_flies_summary = data.frame(
    choice_set = choice_sets
  ) %>% 
    mutate(
      mean_N_flies = as.numeric(mean_N_flies_post)
    ) %>% 
    bind_cols(N_flies_post_percentiles_df)
  
  if(include_all_simulated_responses){
    out = list(
      N_flies_summary = N_flies_summary,
      N_flies_post = N_flies_post
    )
  }else{
    out = list(
      N_flies_summary = N_flies_summary
    )
  }
  
  return(out)
  
}





beta_level_1_draws_hmnl_10samples = create_level_1_betas_hmnl(
  beta_level_0_draws_hmnl, Sigma_level_0_draws_hmnl, n_draws_per_posterior_sample = 10,
  seed = 2023
)





df_to_predict_01 = filtered_exp_10_cs_of_interest_utilities %>% 
  select(folder, choice_set, alternative, pred_ut, R:intensity, no_choice) %>% 
  distinct() %>% 
  select(R:no_choice, folder, choice_set, alternative)


predicted_utilities_df_01 = predict_utilities_design(
  df_to_predict = df_to_predict_01 %>% select(-folder, -choice_set, -alternative),
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = T
)




J = 3
choice_sets = rep(1:(nrow(predicted_utilities_df_01$utilities_summary)/J), each = J)



predicted_N_flies_df_01 = predict_N_flies_design(
  choice_sets = df_to_predict_01$choice_set, 
  utilities_post = predicted_utilities_df_01$utilities_post, 
  n_flies = 80, 
  percentiles = c(5, 10, 25, 50, 75, 90, 95),
  include_all_simulated_responses = T, 
  seed = 2023)

predicted_N_flies_utilities_tibble = predicted_N_flies_df_01$N_flies_summary %>% 
  mutate(alternative = df_to_predict_01$alternative) %>% 
  filter(alternative != "3_none") %>% 
  left_join(utilities_counts_tibble) %>% 
  arrange(group)







pred_counts_exp_10_handpicked_colors_01 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  left_join(
    predicted_N_flies_utilities_tibble %>% 
      select(folder, alternative, group, p10_N_flies, p25_N_flies, p50_N_flies, p75_N_flies, p90_N_flies)
  ) %>% 
  left_join(
    choice_sets_of_interest_catalog %>% 
      select(folder, alternative, color_ix)
  ) %>% 
  mutate(x = paste0(substr(alternative, 1, 1), ") color: ", color_ix)) %>% 
  ggplot(aes(x = x)) +
  geom_jitter(aes(y = count), size = 0.2, width = 0.3, alpha = 0.8) +
  geom_point(aes(y = p50_N_flies), color = "red", size = 1.7) +
  geom_linerange(aes(ymin = p25_N_flies, ymax = p75_N_flies), color = "red", size = 1.3) +
  geom_linerange(aes(ymin = p10_N_flies, ymax = p90_N_flies), color = "red", size = 0.6) +
  # geom_text(aes(y = 1.1*max(max_counts), label = paste0("Median = ", median_counts)),
  #           check_overlap = T, vjust = 1.1) +
  facet_wrap(~group, scales = "free_x") +
  theme_bw() +
  xlab("") +
  ylab("")


ggsave(
  plot = pred_counts_exp_10_handpicked_colors_01,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_01.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)





pred_counts_exp_10_handpicked_colors_02 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  left_join(
    predicted_N_flies_utilities_tibble %>% 
      select(folder, alternative, group, p10_N_flies, p25_N_flies, p50_N_flies, p75_N_flies, p90_N_flies)
  ) %>% 
  left_join(
    choice_sets_of_interest_catalog %>% 
      select(folder, alternative, color_ix)
  ) %>% 
  mutate(x = paste0(substr(alternative, 1, 1), ") color: ", color_ix)) %>% 
  ggplot(aes(x = x)) +
  geom_jitter(aes(y = count), size = 0.2, width = 0.3, alpha = 0.8) +
  geom_point(aes(y = p50_N_flies), color = "red", size = 1.7) +
  geom_linerange(aes(ymin = p25_N_flies, ymax = p75_N_flies), color = "red", size = 1.3) +
  # geom_linerange(aes(ymin = p10_N_flies, ymax = p90_N_flies), color = "red", size = 0.6) +
  # geom_text(aes(y = 1.1*max(max_counts), label = paste0("Median = ", median_counts)),
  #           check_overlap = T, vjust = 1.1) +
  facet_wrap(~group, scales = "free_x") +
  theme_bw() +
  xlab("") +
  ylab("")


ggsave(
  plot = pred_counts_exp_10_handpicked_colors_02,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_02.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)



pred_counts_exp_10_handpicked_colors_03 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  left_join(
    predicted_N_flies_utilities_tibble %>% 
      select(folder, alternative, group, p10_N_flies, p25_N_flies, p50_N_flies, p75_N_flies, p90_N_flies)
  ) %>% 
  left_join(
    choice_sets_of_interest_catalog %>% 
      select(folder, alternative, color_ix)
  ) %>% 
  mutate(x = paste0(substr(alternative, 1, 1), ") color: ", color_ix)) %>% 
  ggplot(aes(x = x)) +
  geom_jitter(aes(y = count), size = 0.2, width = 0.3, alpha = 0.8) +
  geom_point(aes(y = p50_N_flies), color = "red", size = 2) +
  # geom_linerange(aes(ymin = p25_N_flies, ymax = p75_N_flies), color = "red", size = 1.3) +
  # geom_linerange(aes(ymin = p10_N_flies, ymax = p90_N_flies), color = "red", size = 0.6) +
  # geom_text(aes(y = 1.1*max(max_counts), label = paste0("Median = ", median_counts)),
  #           check_overlap = T, vjust = 1.1) +
  facet_wrap(~group, scales = "free_x") +
  theme_bw() +
  xlab("") +
  ylab("")




ggsave(
  plot = pred_counts_exp_10_handpicked_colors_03,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_03.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)




pred_counts_exp_10_handpicked_colors_04 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  left_join(
    predicted_N_flies_utilities_tibble %>% 
      select(folder, alternative, group, p5_N_flies, p10_N_flies, p25_N_flies, p50_N_flies, p75_N_flies, p90_N_flies, p95_N_flies)
  ) %>% 
  left_join(
    choice_sets_of_interest_catalog %>% 
      select(folder, alternative, color_ix)
  ) %>% 
  mutate(x = paste0(substr(alternative, 1, 1), ") color: ", color_ix)) %>% 
  ggplot(aes(x = x)) +
  geom_jitter(aes(y = count), size = 0.2, width = 0.3, alpha = 0.8) +
  geom_point(aes(y = p50_N_flies), color = "red", size = 2) +
  geom_linerange(aes(ymin = p25_N_flies, ymax = p75_N_flies), color = "red", size = 1.3) +
  geom_linerange(aes(ymin = p10_N_flies, ymax = p90_N_flies), color = "red", size = 0.6) +
  geom_linerange(aes(ymin = p5_N_flies, ymax = p95_N_flies), color = "red", size = 0.2) +
  # geom_text(aes(y = 1.1*max(max_counts), label = paste0("Median = ", median_counts)),
  #           check_overlap = T, vjust = 1.1) +
  facet_wrap(~group, scales = "free_x") +
  theme_bw() +
  xlab("") +
  ylab("")




ggsave(
  plot = pred_counts_exp_10_handpicked_colors_04,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_04.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)











predicted_N_flies_df_02 = predict_N_flies_design(
  choice_sets = df_to_predict_01$choice_set, 
  utilities_post = predicted_utilities_df_01$utilities_post, 
  n_flies = 80, 
  percentiles = c(48, 52, seq(5, 95, by = 5)),
  include_all_simulated_responses = F, 
  seed = 2023)




order_hp_colors_plot = rep(NA_character_, 30)
counter = 1
for(i in 1:5){
  for(j in (i+1):6){
    order_hp_colors_plot[counter] = paste(i, j)
    counter = counter + 1
    order_hp_colors_plot[counter] = paste(j, i)
    counter = counter + 1
  }
}
order_hp_colors_plot

order_hp_colors_plot_df = data.frame(
  color_side = order_hp_colors_plot,
  group_order = 1:30
)


predicted_N_flies_utilities_tibble_02 = predicted_N_flies_df_02$N_flies_summary %>%
  mutate(alternative = df_to_predict_01$alternative) %>%
  filter(alternative != "3_none") %>%
  left_join(utilities_counts_tibble %>% 
              select(choice_set, alternative, folder))


predicted_N_flies_utilities_tibble_02_aux = predicted_N_flies_utilities_tibble_02 %>% 
  select(folder, alternative, starts_with("p")) %>% 
  left_join(
    choice_sets_of_interest_catalog %>% 
      select(folder, alternative, color_ix)
  ) %>% 
  group_by(folder) %>% 
  mutate(color_side = paste(color_ix[1], color_ix[2])) %>% 
  ungroup() %>% 
  left_join(order_hp_colors_plot_df) %>% 
  arrange(group_order, alternative) %>% 
  # mutate(group_order = as.integer(factor(color_side))) %>% 
  mutate(group = sprintf("%02d", group_order)) %>% 
  mutate(side = ifelse(substr(alternative, 1, 1) == "1", "L", "R")) %>% 
  mutate(x2 = paste0(side, ": ", color_ix, ifelse(side == "R", "", "   "))) %>% 
  mutate(group2 = paste0("Validation c.s. ", group_order)) %>% 
  mutate(group2 = fct_reorder(group2, group_order))


predicted_N_flies_utilities_tibble_02_aux2 = filtered_exp_10_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  group_by(choice_set, alternative) %>% 
  mutate(median_n_flies = median(count),
         max_n_flies = max(count)) %>% 
  left_join(predicted_N_flies_utilities_tibble_02_aux %>% 
              select(folder, alternative, x2, group2, p90_N_flies)) %>% 
  ungroup()







width_counts_plots_2 = 22
height_counts_plots_2 = 1.4*width_counts_plots_2
linerange_width = 13.5
color_plots = "#505050"





pred_counts_exp_10_handpicked_colors_intervals_04 = predicted_N_flies_utilities_tibble_02_aux2 %>% 
  
  ggplot(aes(x = x2)) +
  
  geom_linerange(
    aes(ymin = p10_N_flies, ymax = p90_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.1,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p20_N_flies, ymax = p80_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux
  ) +
  geom_linerange(
    aes(ymin = p25_N_flies, ymax = p75_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p30_N_flies, ymax = p70_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux
  ) +
  geom_linerange(
    aes(ymin = p45_N_flies, ymax = p55_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.6,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p48_N_flies, ymax = p52_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.9,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  
  
  stat_summary(aes(y = median_n_flies), 
               fun = median, geom = "crossbar", 
               color = "black",
               width = 1, 
               size = 0.3) +
  
  geom_jitter(aes(y = count), size = 0.4, width = 0.3, alpha = 0.8) +
  
  geom_text(aes(x = x2, y = 0.9*plot_height, label = paste0("Median:\n", round(median_n_flies))),
            check_overlap = F, 
            vjust = 1,
            inherit.aes = F,
            size = 3,
            data = predicted_N_flies_utilities_tibble_02_aux2 %>% 
              select(group2, x2, median_n_flies, p90_N_flies) %>% 
              distinct() %>% 
              group_by(group2) %>% 
              mutate(plot_height = max(p90_N_flies)) %>% 
              ungroup()
  ) +
  
  facet_wrap(~group2, scales = "free", ncol = 6) +
  theme_bw() +
  xlab("") +
  ylab("")


ggsave(
  plot = pred_counts_exp_10_handpicked_colors_intervals_04,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_intervals_04.png"),
  width = width_counts_plots_2, height = height_counts_plots_2, units = "cm", dpi = 300
)

knitr::plot_crop(here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_intervals_04.png"))






pred_counts_exp_10_handpicked_colors_intervals_03 = predicted_N_flies_utilities_tibble_02_aux2 %>% 
  
  ggplot(aes(x = x2)) +
  
  geom_linerange(
    aes(ymin = p10_N_flies, ymax = p90_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.1,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p20_N_flies, ymax = p80_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux
  ) +
  geom_linerange(
    aes(ymin = p25_N_flies, ymax = p75_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p30_N_flies, ymax = p70_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux
  ) +
  geom_linerange(
    aes(ymin = p45_N_flies, ymax = p55_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.6,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p48_N_flies, ymax = p52_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.9,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  
  
  stat_summary(aes(y = median_n_flies), 
               fun = median, geom = "crossbar", 
               color = "black",
               width = 1, 
               size = 0.3) +
  
  geom_jitter(aes(y = count), size = 0.4, width = 0.4, alpha = 1) +
  
  facet_wrap(~group2, scales = "free", ncol = 6) +
  theme_bw() +
  xlab("") +
  ylab("")


ggsave(
  plot = pred_counts_exp_10_handpicked_colors_intervals_03,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_intervals_03.png"),
  width = width_counts_plots_2, height = height_counts_plots_2, units = "cm", dpi = 300
)

knitr::plot_crop(here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_intervals_03.png"))







pred_counts_exp_10_handpicked_colors_intervals_02 = predicted_N_flies_utilities_tibble_02_aux2 %>% 
  
  ggplot(aes(x = x2)) +
  
  geom_linerange(
    aes(ymin = p10_N_flies, ymax = p90_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.1,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p20_N_flies, ymax = p80_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux
  ) +
  geom_linerange(
    aes(ymin = p25_N_flies, ymax = p75_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p30_N_flies, ymax = p70_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux
  ) +
  geom_linerange(
    aes(ymin = p45_N_flies, ymax = p55_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.6,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p48_N_flies, ymax = p52_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.9,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  
  
  stat_summary(aes(y = p50_N_flies), 
               fun = median, geom = "crossbar", 
               color = "black",
               width = 1,
               data = predicted_N_flies_utilities_tibble_02_aux) +
  
  geom_jitter(aes(y = count), size = 0.4, width = 0.4, alpha = 1) +
  
  facet_wrap(~group2, scales = "free", ncol = 6) +
  theme_bw() +
  xlab("") +
  ylab("")


ggsave(
  plot = pred_counts_exp_10_handpicked_colors_intervals_02,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_intervals_02.png"),
  width = width_counts_plots_2, height = height_counts_plots_2, units = "cm", dpi = 300
)

knitr::plot_crop(here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_intervals_02.png"))










pred_counts_exp_10_handpicked_colors_intervals_01 = predicted_N_flies_utilities_tibble_02_aux2 %>% 
  
  ggplot(aes(x = x2)) +
  
  geom_linerange(
    aes(ymin = p10_N_flies, ymax = p90_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.1,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p20_N_flies, ymax = p80_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux
  ) +
  geom_linerange(
    aes(ymin = p25_N_flies, ymax = p75_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p30_N_flies, ymax = p70_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_utilities_tibble_02_aux
  ) +
  geom_linerange(
    aes(ymin = p45_N_flies, ymax = p55_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.6,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  geom_linerange(
    aes(ymin = p48_N_flies, ymax = p52_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.9,
    data = predicted_N_flies_utilities_tibble_02_aux) +
  
  geom_jitter(aes(y = count), size = 0.4, width = 0.4, alpha = 1) +
  # geom_point(aes(y = p50_N_flies), color = "red", size = 1.7) +
  
  facet_wrap(~group2, scales = "free", ncol = 6) +
  theme_bw() +
  xlab("") +
  ylab("")


ggsave(
  plot = pred_counts_exp_10_handpicked_colors_intervals_01,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_handpicked_colors_intervals_01.png"),
  width = width_counts_plots_2, height = height_counts_plots_2, units = "cm", dpi = 300
)










####### For the other choice sets


choice_sets_not_of_interest = unique(design_confirmatory_ordered$cs_new[61:nrow(design_confirmatory_ordered)])


choice_sets_to_analyze_not_interest_folder = design_confirmatory_with_mixture[choice_sets_not_of_interest, "folder"]
# choice_sets_to_analyze = design_confirmatory_with_mixture[1:28, "folder"]


filtered_exp_10_cs_of_not_interest = counts_japanese %>% 
  filter(experiment == 10, folder %in% choice_sets_to_analyze_not_interest_folder) %>% 
  select(R:no_choice)


df_to_predict_cs_not_of_interest = counts_japanese %>% 
  filter(experiment == 10, folder %in% choice_sets_to_analyze_not_interest_folder) %>% 
  select(folder, choice_set, alternative, R:intensity, no_choice) %>% 
  distinct() %>% 
  select(R:no_choice, folder, choice_set, alternative)


predicted_utilities_df_cs_not_of_interest = predict_utilities_design(
  df_to_predict = df_to_predict_cs_not_of_interest %>% select(-folder, -choice_set, -alternative),
  beta_level_1_draws_hmnl = beta_level_1_draws_hmnl_10samples,
  include_all_utilities = T
)



predicted_N_flies_df_cs_not_of_interest = predict_N_flies_design(
  choice_sets = df_to_predict_cs_not_of_interest$choice_set, 
  utilities_post = predicted_utilities_df_cs_not_of_interest$utilities_post, 
  n_flies = 80, 
  percentiles = c(48, 52, seq(5, 95, by = 5)),
  include_all_simulated_responses = F, 
  seed = 2023)





predicted_N_flies_df_cs_not_of_interest_part1 = predicted_N_flies_df_cs_not_of_interest$N_flies_summary %>%
  mutate(alternative = df_to_predict_cs_not_of_interest$alternative) %>%
  filter(alternative != "3_none") %>% 
  left_join(
    df_to_predict_cs_not_of_interest %>% 
      select(choice_set, folder, alternative)
  ) %>% 
  filter(folder %in% choice_sets_to_analyze_not_interest_folder[1:30])



predicted_N_flies_df_cs_not_of_interest_part2 = predicted_N_flies_df_cs_not_of_interest$N_flies_summary %>%
  mutate(alternative = df_to_predict_cs_not_of_interest$alternative) %>%
  filter(alternative != "3_none") %>% 
  left_join(
    df_to_predict_cs_not_of_interest %>% 
      select(choice_set, folder, alternative)
  ) %>% 
  filter(folder %in% choice_sets_to_analyze_not_interest_folder[31:58])




predicted_N_flies_df_cs_not_of_interest_aux2_part1 = filtered_exp_10_cs_of_not_interest %>% 
  filter(alternative != "3_none") %>% 
  group_by(choice_set, alternative) %>% 
  mutate(median_n_flies = median(count),
         max_n_flies = max(count)) %>% 
  filter(folder %in% choice_sets_to_analyze_not_interest_folder[1:30])



predicted_N_flies_df_cs_not_of_interest_aux2_part2 = filtered_exp_10_cs_of_not_interest %>% 
  filter(alternative != "3_none") %>% 
  group_by(choice_set, alternative) %>% 
  mutate(median_n_flies = median(count),
         max_n_flies = max(count)) %>% 
  filter(folder %in% choice_sets_to_analyze_not_interest_folder[31:58])





### Plot of part 1

pred_counts_exp_10_cs_not_interest_intervals_part1_01 = predicted_N_flies_df_cs_not_of_interest_aux2_part1 %>% 
  
  ggplot(aes(x = alternative)) +
  
  geom_linerange(
    aes(ymin = p10_N_flies, ymax = p90_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.1,
    data = predicted_N_flies_df_cs_not_of_interest_part1) +
  geom_linerange(
    aes(ymin = p20_N_flies, ymax = p80_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_df_cs_not_of_interest_part1
  ) +
  geom_linerange(
    aes(ymin = p25_N_flies, ymax = p75_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_df_cs_not_of_interest_part1) +
  geom_linerange(
    aes(ymin = p30_N_flies, ymax = p70_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_df_cs_not_of_interest_part1
  ) +
  geom_linerange(
    aes(ymin = p45_N_flies, ymax = p55_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.6,
    data = predicted_N_flies_df_cs_not_of_interest_part1) +
  geom_linerange(
    aes(ymin = p48_N_flies, ymax = p52_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.9,
    data = predicted_N_flies_df_cs_not_of_interest_part1) +
  
  
  stat_summary(aes(y = median_n_flies), 
               fun = median, geom = "crossbar", 
               color = "black",
               width = 1, 
               size = 0.3) +
  
  geom_jitter(aes(y = count), size = 0.4, width = 0.3, alpha = 0.8) +
  
  # geom_text(aes(x = alternative, y = 0.9*plot_height, label = paste0("Median:\n", round(median_n_flies))),
  #           check_overlap = F,
  #           vjust = 1,
  #           inherit.aes = F,
  #           size = 3,
  #           data = predicted_N_flies_utilities_tibble_02_aux2 %>%
  #             select(group2, x2, median_n_flies, p90_N_flies) %>%
  #             distinct() %>%
  #             group_by(group2) %>%
  #             mutate(plot_height = max(p90_N_flies)) %>%
  #             ungroup()
  # ) +

facet_wrap(~choice_set, scales = "free", ncol = 6) +
  theme_bw() +
  xlab("") +
  ylab("")




ggsave(
  plot = pred_counts_exp_10_cs_not_interest_intervals_part1_01,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_cs_not_interest_intervals_part1_01.png"),
  width = width_counts_plots_2, height = height_counts_plots_2, units = "cm", dpi = 300
)




### Plot of part 2

pred_counts_exp_10_cs_not_interest_intervals_part2_01 = predicted_N_flies_df_cs_not_of_interest_aux2_part2 %>% 
  
  ggplot(aes(x = alternative)) +
  
  geom_linerange(
    aes(ymin = p10_N_flies, ymax = p90_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.1,
    data = predicted_N_flies_df_cs_not_of_interest_part2) +
  geom_linerange(
    aes(ymin = p20_N_flies, ymax = p80_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_df_cs_not_of_interest_part2
  ) +
  geom_linerange(
    aes(ymin = p25_N_flies, ymax = p75_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_df_cs_not_of_interest_part2) +
  geom_linerange(
    aes(ymin = p30_N_flies, ymax = p70_N_flies), 
    color = color_plots, size = linerange_width, alpha = 0.2,
    data = predicted_N_flies_df_cs_not_of_interest_part2
  ) +
  geom_linerange(
    aes(ymin = p45_N_flies, ymax = p55_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.6,
    data = predicted_N_flies_df_cs_not_of_interest_part2) +
  geom_linerange(
    aes(ymin = p48_N_flies, ymax = p52_N_flies),
    color = color_plots, size = linerange_width, alpha = 0.9,
    data = predicted_N_flies_df_cs_not_of_interest_part2) +
  
  
  stat_summary(aes(y = median_n_flies), 
               fun = median, geom = "crossbar", 
               color = "black",
               width = 1, 
               size = 0.3) +
  
  geom_jitter(aes(y = count), size = 0.4, width = 0.3, alpha = 0.8) +
  
  # geom_text(aes(x = alternative, y = 0.9*plot_height, label = paste0("Median:\n", round(median_n_flies))),
  #           check_overlap = F,
  #           vjust = 1,
  #           inherit.aes = F,
  #           size = 3,
  #           data = predicted_N_flies_utilities_tibble_02_aux2 %>%
  #             select(group2, x2, median_n_flies, p90_N_flies) %>%
  #             distinct() %>%
  #             group_by(group2) %>%
  #             mutate(plot_height = max(p90_N_flies)) %>%
#             ungroup()
# ) +

facet_wrap(~choice_set, scales = "free", ncol = 6) +
  theme_bw() +
  xlab("") +
  ylab("")


ggsave(
  plot = pred_counts_exp_10_cs_not_interest_intervals_part2_01,
  filename = here("japanese_fly/out/plots/pred_counts_exp_10_cs_not_interest_intervals_part2_01.png"),
  width = width_counts_plots_2, height = height_counts_plots_2, units = "cm", dpi = 300
)



# Predictions on which alternative is more attractive according to the median number of flies
# and the predicted number of flies

predicted_N_flies_utilities_tibble_02_aux2 %>% 
  select(folder, choice_set, alternative, group2, median_n_flies) %>% 
  distinct() %>% 
  left_join(
    predicted_N_flies_utilities_tibble_02
  ) %>% 
  group_by(folder, choice_set, group2) %>% 
  mutate(
    n_left_real = first(median_n_flies),
    n_right_real = last(median_n_flies),
    n_left_pred_p45 = first(p45_N_flies),
    n_right_pred_p45 = last(p45_N_flies),
    n_left_pred_p55 = first(p55_N_flies),
    n_right_pred_p55 = last(p55_N_flies)
  ) %>% 
  select(starts_with("n_")) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(pred_correct = ifelse(
    (n_left_real > n_right_real & n_left_pred_p45 > n_right_pred_p45) |
      (n_left_real > n_right_real & n_left_pred_p55 > n_right_pred_p55) | 
      (n_left_real < n_right_real & n_left_pred_p45 < n_right_pred_p45) |
      (n_left_real < n_right_real & n_left_pred_p55 < n_right_pred_p55) | 
      (n_left_real == n_right_real & n_left_pred_p55 == n_right_pred_p55) |
      (n_left_real == n_right_real & n_left_pred_p45 == n_right_pred_p45), 
    1, 0)) %>% 
  pull(pred_correct) %>% 
  sum()

# 26 cs out of 30




predicted_N_flies_df_cs_not_of_interest_aux2_part1 %>% 
  select(folder, choice_set, alternative, median_n_flies) %>% 
  distinct() %>% 
  left_join(
    predicted_N_flies_df_cs_not_of_interest_part1
  ) %>% 
  group_by(folder, choice_set) %>% 
  mutate(
    n_left_real = first(median_n_flies),
    n_right_real = last(median_n_flies),
    n_left_pred_p45 = first(p45_N_flies),
    n_right_pred_p45 = last(p45_N_flies),
    n_left_pred_p55 = first(p55_N_flies),
    n_right_pred_p55 = last(p55_N_flies)
  ) %>% 
  select(starts_with("n_")) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(pred_correct = ifelse(
    (n_left_real > n_right_real & n_left_pred_p45 > n_right_pred_p45) |
      (n_left_real > n_right_real & n_left_pred_p55 > n_right_pred_p55) | 
      (n_left_real < n_right_real & n_left_pred_p45 < n_right_pred_p45) |
      (n_left_real < n_right_real & n_left_pred_p55 < n_right_pred_p55) | 
      (n_left_real == n_right_real & n_left_pred_p55 == n_right_pred_p55) |
      (n_left_real == n_right_real & n_left_pred_p45 == n_right_pred_p45), 
    1, 0)) %>% 
  pull(pred_correct) %>% 
  sum()

# 27 cs out of 30







predicted_N_flies_df_cs_not_of_interest_aux2_part2 %>% 
  select(folder, choice_set, alternative, median_n_flies) %>% 
  distinct() %>% 
  left_join(
    predicted_N_flies_df_cs_not_of_interest_part2
  ) %>% 
  group_by(folder, choice_set) %>% 
  mutate(
    n_left_real = first(median_n_flies),
    n_right_real = last(median_n_flies),
    n_left_pred_p45 = first(p45_N_flies),
    n_right_pred_p45 = last(p45_N_flies),
    n_left_pred_p55 = first(p55_N_flies),
    n_right_pred_p55 = last(p55_N_flies)
  ) %>% 
  select(starts_with("n_")) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(pred_correct = ifelse(
    (n_left_real > n_right_real & n_left_pred_p45 > n_right_pred_p45) |
      (n_left_real > n_right_real & n_left_pred_p55 > n_right_pred_p55) | 
      (n_left_real < n_right_real & n_left_pred_p45 < n_right_pred_p45) |
      (n_left_real < n_right_real & n_left_pred_p55 < n_right_pred_p55) | 
      (n_left_real == n_right_real & n_left_pred_p55 == n_right_pred_p55) |
      (n_left_real == n_right_real & n_left_pred_p45 == n_right_pred_p45), 
    1, 0)) %>% 
  pull(pred_correct) %>% 
  sum()

# 20 cs out of 28



# 26 + 27 + 20 = 73 out of 88 choice sets
# 73/88 = 0.818181
# In 82% of choice sets our predictions matched which would be more attractive

# The numbers are almost identical using 48 and 52 percentiles instead of 45 and 55:
# 25 + 27 + 21 = 73 out of 88 choice sets





##### 100% UV


folders_with_100uv = counts_japanese %>% 
  filter(UV == 1, no_choice == 0) %>% 
  pull(folder) %>% 
  unique()


counts_japanese %>% 
  filter(folder %in% folders_with_100uv)


counts_japanese %>% 
  filter(folder %in% folders_with_100uv) %>% 
  filter(UV == 1) %>% 
  group_by(experiment, folder, alternative) %>% 
  summarize(median_n_flies = median(count),
            max_n_flies = max(count))



uv_100_tibble_to_predict = tibble(
  R = 0, O = 0, Y = 0, G = 0, B = 0, P = 0, UV = 1, intensity = seq(-1, 1, by = 0.05)
) %>% 
  mutate(no_choice = 0)


uv_100_predicted_utilities_samples = predict_utilities_design(
  df_to_predict = uv_100_tibble_to_predict,
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = T,
  n_draws_per_posterior_sample = 10,
  seed = 2023
)



uv_100_predicted_utilities_samples$utilities_summary %>% 
  ggplot(aes(x = intensity, y = p50, ymin = p10, ymax = p90)) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  theme_bw() +
  ylab("Predicted utility")



best_color_diff_intensities_tibble = exp_10_handpicked_colors_indices_and_utilities[1,] %>% 
  select(R:UV) %>% 
  bind_cols(intensity = seq(-1, 1, by = 0.05)) %>% 
  mutate(no_choice = 0)



best_color_diff_intensities_predicted_utilities_samples = predict_utilities_design(
  df_to_predict = best_color_diff_intensities_tibble,
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = F,
  n_draws_per_posterior_sample = 10,
  seed = 2023
)



best_color_diff_intensities_predicted_utilities_samples$utilities_summary %>% 
  ggplot(aes(x = intensity, y = p50, ymin = p10, ymax = p90)) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  theme_bw() +
  ylab("Predicted utility")




utility_intensity_plots_width = 20
utility_intensity_plots_height = 10

best_color_100uv_intensity_plot = best_color_diff_intensities_predicted_utilities_samples$utilities_summary %>% 
  mutate(Color = "Best predicted color") %>% 
  bind_rows(
    uv_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% UV")
  ) %>% 
  ggplot(aes(x = intensity, y = p50, ymin = p10, ymax = p90, color = Color, fill = Color)) +
  geom_line(size = 1.2) +
  geom_ribbon(alpha = 0.2, linetype = 0) +
  theme_bw() +
  xlab("Intensity") +
  ylab("Predicted utility") +
  scale_fill_manual(values = c("#5f4b8b", "#28bd5f")) +
  scale_color_manual(values = c("#5f4b8b", "#28bd5f"))
  

ggsave(
  plot = best_color_100uv_intensity_plot,
  filename = here("japanese_fly/out/plots/best_color_100uv_intensity_plot_01.png"),
  width = utility_intensity_plots_width, height = utility_intensity_plots_height, 
  units = "cm", dpi = 300
)










r_100_predicted_utilities_samples = predict_utilities_design(
  df_to_predict = tibble(
    R = 1, O = 0, Y = 0, G = 0, B = 0, P = 0, UV = 0, intensity = seq(-1, 1, by = 0.05)
  ) %>% 
    mutate(no_choice = 0),
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = F,
  n_draws_per_posterior_sample = 10,
  seed = 2023
)


o_100_predicted_utilities_samples = predict_utilities_design(
  df_to_predict = tibble(
    R = 0, O = 1, Y = 0, G = 0, B = 0, P = 0, UV = 0, intensity = seq(-1, 1, by = 0.05)
  ) %>% 
    mutate(no_choice = 0),
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = F,
  n_draws_per_posterior_sample = 10,
  seed = 2023
)


y_100_predicted_utilities_samples = predict_utilities_design(
  df_to_predict = tibble(
    R = 0, O = 0, Y = 1, G = 0, B = 0, P = 0, UV = 0, intensity = seq(-1, 1, by = 0.05)
  ) %>% 
    mutate(no_choice = 0),
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = F,
  n_draws_per_posterior_sample = 10,
  seed = 2023
)


g_100_predicted_utilities_samples = predict_utilities_design(
  df_to_predict = tibble(
    R = 0, O = 0, Y = 0, G = 1, B = 0, P = 0, UV = 0, intensity = seq(-1, 1, by = 0.05)
  ) %>% 
    mutate(no_choice = 0),
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = F,
  n_draws_per_posterior_sample = 10,
  seed = 2023
)


b_100_predicted_utilities_samples = predict_utilities_design(
  df_to_predict = tibble(
    R = 0, O = 0, Y = 0, G = 0, B = 1, P = 0, UV = 0, intensity = seq(-1, 1, by = 0.05)
  ) %>% 
    mutate(no_choice = 0),
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = F,
  n_draws_per_posterior_sample = 10,
  seed = 2023
)


p_100_predicted_utilities_samples = predict_utilities_design(
  df_to_predict = tibble(
    R = 0, O = 0, Y = 0, G = 0, B = 0, P = 1, UV = 0, intensity = seq(-1, 1, by = 0.05)
  ) %>% 
    mutate(no_choice = 0),
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = F,
  n_draws_per_posterior_sample = 10,
  seed = 2023
)










colors_plot_100percent = c(
  "blue",
  "dark green",
  "orange",
  "purple",
  "red",
  "#5f4b8b",
  "yellow",
  "#28bd5f"
)


best_color_100other_colors_intensity_plot_01 = b_100_predicted_utilities_samples$utilities_summary %>% 
  mutate(Color = "100% blue") %>% 
  bind_rows(
    g_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% green")
  ) %>% 
  bind_rows(
    o_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% orange")
  ) %>% 
  bind_rows(
    p_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% purple")
  ) %>% 
  bind_rows(
    r_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% red")
  ) %>% 
  bind_rows(
    uv_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% UV")
  ) %>% 
  bind_rows(
    y_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% yellow")
  ) %>% 
  bind_rows(
    best_color_diff_intensities_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "Best predicted color")
  ) %>% 
  # mutate(ix = 1:nrow(.)) %>% 
  # mutate(fct_reorder(Color, ix)) %>% 
  ggplot(aes(x = intensity, y = p50, ymin = p10, ymax = p90, color = Color, fill = Color)) +
  geom_line(size = 1.2) +
  geom_ribbon(alpha = 0.2, linetype = 0) +
  theme_bw() +
  xlab("Intensity") +
  ylab("Predicted utility") +
  scale_color_manual(values = colors_plot_100percent) +
  scale_fill_manual(values = colors_plot_100percent)

ggsave(
  plot = best_color_100other_colors_intensity_plot_01,
  filename = here("japanese_fly/out/plots/best_color_100other_colors_intensity_plot_01.png"),
  width = utility_intensity_plots_width, height = utility_intensity_plots_height, 
  units = "cm", dpi = 300
)



best_color_100other_colors_intensity_plot_02 = b_100_predicted_utilities_samples$utilities_summary %>% 
  mutate(Color = "100% blue") %>% 
  bind_rows(
    g_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% green")
  ) %>% 
  bind_rows(
    o_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% orange")
  ) %>% 
  bind_rows(
    p_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% purple")
  ) %>% 
  bind_rows(
    r_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% red")
  ) %>% 
  bind_rows(
    uv_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% UV")
  ) %>% 
  bind_rows(
    y_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% yellow")
  ) %>% 
  bind_rows(
    best_color_diff_intensities_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "Best predicted color")
  ) %>% 
  # mutate(ix = 1:nrow(.)) %>% 
  # mutate(fct_reorder(Color, ix)) %>% 
  ggplot(aes(x = intensity, y = p50, ymin = p50 - sd, ymax = p50 + sd, color = Color, fill = Color)) +
  geom_line(size = 1.2) +
  geom_ribbon(alpha = 0.2, linetype = 0) +
  theme_bw() +
  xlab("Intensity") +
  ylab("Predicted utility") +
  scale_color_manual(values = colors_plot_100percent) +
  scale_fill_manual(values = colors_plot_100percent)


ggsave(
  plot = best_color_100other_colors_intensity_plot_02,
  filename = here("japanese_fly/out/plots/best_color_100other_colors_intensity_plot_02.png"),
  width = utility_intensity_plots_width, height = utility_intensity_plots_height, 
  units = "cm", dpi = 300
)



best_color_100other_colors_intensity_plot_03 = b_100_predicted_utilities_samples$utilities_summary %>% 
  mutate(Color = "100% blue") %>% 
  bind_rows(
    g_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% green")
  ) %>% 
  bind_rows(
    o_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% orange")
  ) %>% 
  bind_rows(
    p_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% purple")
  ) %>% 
  bind_rows(
    r_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% red")
  ) %>% 
  bind_rows(
    uv_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% UV")
  ) %>% 
  bind_rows(
    y_100_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "100% yellow")
  ) %>% 
  bind_rows(
    best_color_diff_intensities_predicted_utilities_samples$utilities_summary %>% 
      mutate(Color = "Best predicted color")
  ) %>% 
  # mutate(ix = 1:nrow(.)) %>% 
  # mutate(fct_reorder(Color, ix)) %>% 
  ggplot(aes(x = intensity, y = p50, ymin = p50 - sd, ymax = p50 + sd, color = Color, fill = Color)) +
  geom_line(size = 1.2) +
  # geom_ribbon(alpha = 0.2, linetype = 0) +
  theme_bw() +
  xlab("Intensity") +
  ylab("Predicted utility") +
  scale_color_manual(values = colors_plot_100percent) +
  scale_fill_manual(values = colors_plot_100percent)


ggsave(
  plot = best_color_100other_colors_intensity_plot_03,
  filename = here("japanese_fly/out/plots/best_color_100other_colors_intensity_plot_03.png"),
  width = utility_intensity_plots_width, height = utility_intensity_plots_height, 
  units = "cm", dpi = 300
)






# Distribution of utilities in each experiment

choice_sets_all_exps = counts_japanese %>% 
  filter(no_choice != 1) %>% 
  select(experiment, choice_set, folder, alternative, R:no_choice) %>% 
  distinct()

betas_level_0_summary

utilities_choice_sets_all_exps = (choice_sets_all_exps %>% 
                                    select(R:no_choice) %>% 
                                    create_model_matrix_second_order_scheffe())$X %*% betas_level_0_summary$mean

choice_sets_all_exps$exp_utility = utilities_choice_sets_all_exps

choice_sets_all_exps %>% 
  ggplot() +
  geom_boxplot(aes(x = experiment, y = exp_utility, group = experiment)) +
  xlab("Experiment") +
  ylab("Expected utility") +
  theme_bw()










#### Attempts at some other things



exp_10_diff_counts = counts_japanese %>% 
  filter(experiment == 10, alternative != "3_none") %>% 
  select(folder, image, alternative, count) %>% 
  group_by(folder, image) %>% 
  mutate(count_right = lead(count)) %>% 
  ungroup() %>% 
  filter(complete.cases(.)) %>% 
  mutate(right_minus_left = count_right - count) %>% 
  select(folder, image, count_left = count, count_right, right_minus_left)



exp_10_colors_unique_to_predict = counts_japanese %>% 
  filter(experiment == 10, alternative != "3_none") %>% 
  select(folder, alternative, R:no_choice) %>% 
  distinct()












exp_10_colors_unique_to_predict_list <- create_model_matrix_second_order_scheffe(exp_10_colors_unique_to_predict %>% select(-folder, -alternative))


exp_10_colors_unique_predicted_point_utilities = exp_10_colors_unique_to_predict_list$X %*% betas_level_0_summary$mean


exp_10_colors_unique_to_predict_utilities = exp_10_colors_unique_to_predict %>% 
  mutate(pred_ut = as.numeric(exp_10_colors_unique_predicted_point_utilities))


exp10_predicted_point_utilities_per_cs = exp_10_colors_unique_to_predict_utilities %>% 
  group_by(folder) %>% 
  mutate(left_ut = exp(first(pred_ut)),
         right_ut = exp(last(pred_ut))) %>% 
  select(folder, left_ut, right_ut) %>% 
  distinct()





exp_10_diff_counts_utilities_preds = exp_10_diff_counts %>% 
  left_join(exp10_predicted_point_utilities_per_cs) %>% 
  mutate(model_right = ifelse(right_ut > left_ut & right_minus_left > 0, 1, 0))

exp_10_diff_counts_utilities_preds %>% 
  group_by(folder) %>% 
  summarize(
    prop_model_right = mean(model_right),
    median_n_flies = median(count_left + count_right),
    max_n_flies = max(count_left + count_right)
  ) %>% 
  arrange(desc(prop_model_right)) %>% 
  left_join(exp10_predicted_point_utilities_per_cs) %>% 
  mutate(dif_ut = right_ut - left_ut)











# Compute probabilities of right > left

exp_10_colors_unique_predicted_utilities_samples = predict_utilities_design(
  df_to_predict = exp_10_colors_unique_to_predict %>% select(-folder, -alternative),
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = T,
  n_draws_per_posterior_sample = 10,
  seed = 2023
)


exp_10_colors_unique_predicted_utilities_samples$utilities_post %>% dim()

exp_10_right_bigger_than_left_samples = matrix(
  NA_integer_,
  nrow = nrow(exp_10_colors_unique_predicted_utilities_samples$utilities_post)/2,
  ncol = ncol(exp_10_colors_unique_predicted_utilities_samples$utilities_post),
)


for(r in 1:ncol(exp_10_right_bigger_than_left_samples)){
  if(r %% 5000 == 1) cat("Sample", r, "\n")
  for(j in 1:nrow(exp_10_right_bigger_than_left_samples)){
    left_alt = 2*j - 1
    right_alt = 2*j
    res = exp_10_colors_unique_predicted_utilities_samples$utilities_post[right_alt, r] >
      exp_10_colors_unique_predicted_utilities_samples$utilities_post[left_alt, r]
    
    exp_10_right_bigger_than_left_samples[j, r] = as.integer(res)
  }
}


exp_10_probs_right_better_left = data.frame(
  folder = unique(exp_10_colors_unique_to_predict$folder),
  prob_right_better_than_left = apply(exp_10_right_bigger_than_left_samples, 1, mean)
)



# exp_10_diff_counts_utilities_preds %>% 
#   group_by(folder) %>% 
#   summarize(
#     prop_model_right = mean(model_right),
#     median_n_flies = median(count_left + count_right),
#     max_n_flies = max(count_left + count_right)
#   ) %>% 
#   arrange(desc(prop_model_right)) %>% 
#   left_join(exp10_predicted_point_utilities_per_cs) %>% 
#   mutate(dif_ut = right_ut - left_ut) %>% 
#   left_join(exp_10_probs_right_better_left) %>% 
#   View()


exp_10_diff_counts %>% 
  group_by(folder) %>% 
  summarize(
    # prop_model_right = mean(model_right),
    median_n_flies = median(count_left + count_right),
    max_n_flies = max(count_left + count_right)
  ) %>% 
  # arrange(desc(prop_model_right)) %>% 
  # left_join(exp10_predicted_point_utilities_per_cs) %>% 
  # mutate(dif_ut = right_ut - left_ut) %>% 
  left_join(exp_10_probs_right_better_left) %>% 
  View()






exp_10_diff_proportions_probs = exp_10_diff_counts %>% 
  group_by(folder) %>% 
  summarize(
    median_left = median(count_left),
    median_right = median(count_right),
    median_n_total_flies = median(count_left + count_right),
    prop_right_more_than_left = mean(count_right > count_left)
  ) %>% 
  left_join(exp_10_probs_right_better_left) %>% 
  mutate(pred_minus_real = prob_right_better_than_left - prop_right_more_than_left)

exp_10_diff_proportions_probs

exp_10_diff_proportions_probs %>% 
  arrange(desc(prop_right_more_than_left)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity") +
  geom_point(aes(y = prob_right_better_than_left))



exp_10_diff_proportions_probs %>% 
  arrange(desc(prob_right_better_than_left)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity") +
  geom_point(aes(y = prob_right_better_than_left))




exp_10_diff_proportions_probs %>% 
  arrange(desc(median_n_total_flies)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  # geom_bar(aes(y = prop_right_more_than_left), stat = "identity") +
  geom_point(aes(y = pred_minus_real, size = median_n_total_flies))




exp_10_diff_proportions_probs %>% 
  arrange(desc(prop_right_more_than_left)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  # geom_bar(aes(y = prop_right_more_than_left), stat = "identity") +
  geom_point(aes(y = pred_minus_real, size = median_n_total_flies))

















# Compute probabilities of N_right > N_left

n_images_per_cs = 120

exp_10_df_to_predict = counts_japanese %>% 
  filter(experiment == 10) %>% 
  select(folder, alternative, R:no_choice) %>% 
  distinct()

exp_10_df_to_predict_utilities_samples = predict_utilities_design(
  df_to_predict = exp_10_df_to_predict %>% select(-folder, -alternative),
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = T,
  n_draws_per_posterior_sample = n_images_per_cs,
  seed = 2023
)


choice_sets_exp_10 = rep(1:(nrow(exp_10_df_to_predict_utilities_samples$utilities_summary)/3), each = 3)

predicted_N_flies_exp_10 = predict_N_flies_design(
  choice_sets = choice_sets_exp_10,
  utilities_post = exp_10_df_to_predict_utilities_samples$utilities_post, 
  n_flies = 80, 
  percentiles = c(5, 10, 25, 50, 75, 90, 95),
  include_all_simulated_responses = T, 
  seed = 2023)

predicted_N_flies_exp_10$N_flies_summary %>% 
  mutate(alternative = exp_10_df_to_predict_utilities_samples$alternative)



exp_10_count_right_more_than_left_samples = array(
  NA_integer_,
  c(dim(predicted_N_flies_exp_10$N_flies_post)[1]/n_images_per_cs, 
    n_images_per_cs,
    dim(predicted_N_flies_exp_10$N_flies_post)[3]
  )
  
)


r2 = 0
for(i in 1:dim(exp_10_count_right_more_than_left_samples)[2]){
  for(r in 1:dim(exp_10_count_right_more_than_left_samples)[1]){
    r2 = r2 + 1
    if(r2 %% 5000 == 1) cat("Sample", r2, "\n")
    for(s in 1:dim(exp_10_count_right_more_than_left_samples)[3]){
      
      left_alt = 1
      right_alt = 2
      res = predicted_N_flies_exp_10$N_flies_post[r2, right_alt, s] >
        predicted_N_flies_exp_10$N_flies_post[r2, left_alt, s]
      
      exp_10_count_right_more_than_left_samples[r, i, s] = as.integer(res)
    }
  }
}



# Distribution of proportion of images where right alt has more flies than left alt
exp_10_distr_proportion_of_images_r_bigger_l = apply(exp_10_count_right_more_than_left_samples, c(1, 3), mean)

# Summary
exp_10_distr_proportion_of_images_r_bigger_l_summary = t(apply(exp_10_distr_proportion_of_images_r_bigger_l, 2, quantile, c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))) %>% 
  as.data.frame() %>% 
  set_names(c("p025", "p05", "p10", "p50", "p90", "p95", "p975")) %>% 
  mutate(folder = unique(exp_10_df_to_predict$folder))

head(exp_10_distr_proportion_of_images_r_bigger_l_summary)









exp_10_diff_proportions_probs_percentiles = exp_10_diff_counts %>% 
  group_by(folder) %>% 
  summarize(
    median_left = median(count_left),
    median_right = median(count_right),
    median_n_total_flies = median(count_left + count_right),
    prop_right_more_than_left = mean(count_right > count_left)
  ) %>% 
  left_join(exp_10_distr_proportion_of_images_r_bigger_l_summary)

exp_10_diff_proportions_probs_percentiles

exp_10_diff_proportions_probs_percentiles %>% 
  arrange(desc(prop_right_more_than_left)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity") +
  geom_point(aes(y = p50)) +
  geom_linerange(aes(ymin = p025, ymax = p975))





exp_10_diff_proportions_probs_percentiles %>% 
  arrange(desc(p50)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity") +
  geom_point(aes(y = p50)) +
  geom_linerange(aes(ymin = p025, ymax = p975))

































# Distribution of probs of utilities

n_images_per_cs = 120

exp_10_repeated_utilities_samples = predict_utilities_design(
  df_to_predict = exp_10_colors_unique_to_predict %>% select(-folder, -alternative),
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl,
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = T,
  n_draws_per_posterior_sample = 2*n_images_per_cs,
  seed = 2023
)




exp_10_right_better_than_left_repeated_samples = array(
  NA_integer_,
  c(ncol(exp_10_repeated_utilities_samples$utilities_post)/n_images_per_cs, 
    n_images_per_cs,
    nrow(exp_10_repeated_utilities_samples$utilities_post)/2
  )
  
)


r2 = 0
for(i in 1:dim(exp_10_right_better_than_left_repeated_samples)[2]){
  for(r in 1:dim(exp_10_right_better_than_left_repeated_samples)[1]){
    r2 = r2 + 1
    if(r2 %% 5000 == 1) cat("Sample", r2, "\n")
    for(s in 1:dim(exp_10_right_better_than_left_repeated_samples)[3]){
      
      left_alt = 2*s-1
      right_alt = 2*s
      res = exp_10_repeated_utilities_samples$utilities_post[right_alt, r2] >
        exp_10_repeated_utilities_samples$utilities_post[left_alt, r2]
      
      exp_10_right_better_than_left_repeated_samples[r, i, s] = as.integer(res)
    }
  }
}

# exp_10_repeated_utilities_samples$utilities_post <- NULL


# Distribution of proportion of images where right alt has a higher utility than left alt
exp_10_distr_proportion_of_images_r_better_l = apply(exp_10_right_better_than_left_repeated_samples, c(1, 3), mean)

# Summary
exp_10_distr_proportion_of_images_r_better_l_summary = t(apply(exp_10_distr_proportion_of_images_r_better_l, 2, quantile, c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))) %>% 
  as.data.frame() %>% 
  set_names(c("p025", "p05", "p10", "p50", "p90", "p95", "p975")) %>% 
  mutate(folder = unique(exp_10_colors_unique_to_predict$folder),
         ix = 1:nrow(.))

head(exp_10_distr_proportion_of_images_r_better_l_summary)







exp_10_repeated_utilities_samples_left_summary = exp_10_repeated_utilities_samples$utilities_summary %>% 
  slice(seq(1, nrow(.), by = 2)) %>% 
  mutate(left_ut_p10 = exp(p10), left_ut_p50 = exp(p50), left_ut_p90 = exp(p90)) %>% 
  select(starts_with("left"))


exp_10_repeated_utilities_samples_right_summary = exp_10_repeated_utilities_samples$utilities_summary %>% 
  slice(seq(2, nrow(.), by = 2)) %>% 
  mutate(right_ut_p10 = exp(p10), right_ut_p50 = exp(p50), right_ut_p90 = exp(p90)) %>% 
  select(starts_with("right"))


exp_10_diff_proportions_probs_percentiles_utilities = exp_10_diff_counts %>% 
  group_by(folder) %>% 
  summarize(
    median_left = median(count_left),
    median_right = median(count_right),
    median_n_total_flies = median(count_left + count_right),
    max_n_total_flies = max(count_left + count_right),
    median_rel_diff = median(abs(count_right - count_left)/(count_left + count_right), na.rm = T),
    prop_right_more_than_left = mean(count_right > count_left)
  ) %>% 
  left_join(exp_10_distr_proportion_of_images_r_better_l_summary) %>% 
  arrange(ix) %>% 
  bind_cols(exp_10_repeated_utilities_samples_left_summary) %>% 
  bind_cols(exp_10_repeated_utilities_samples_right_summary)

exp_10_diff_proportions_probs_percentiles_utilities

exp_10_diff_proportions_probs_percentiles_utilities %>% 
  arrange(desc(prop_right_more_than_left)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity", fill = "dark gray") +
  geom_point(aes(y = p50)) +
  geom_linerange(aes(ymin = p025, ymax = p975)) +
  theme_bw()





exp_10_diff_proportions_probs_percentiles_utilities %>% 
  arrange(desc(p50)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity", fill = "dark gray") +
  geom_point(aes(y = p50)) +
  geom_linerange(aes(ymin = p025, ymax = p975)) +
  theme_bw()





exp_10_diff_proportions_probs_percentiles_utilities %>% 
  mutate(pred_minus_real = abs(p50 - prop_right_more_than_left)) %>% 
  arrange(pred_minus_real) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity", fill = "dark gray") +
  geom_point(aes(y = p50)) +
  geom_linerange(aes(ymin = p025, ymax = p975)) +
  theme_bw()



exp_10_diff_proportions_probs_percentiles_utilities %>% 
  arrange(desc(median_n_total_flies)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity", fill = "dark gray") +
  geom_point(aes(y = p50, size = median_n_total_flies)) +
  geom_linerange(aes(ymin = p025, ymax = p975)) +
  theme_bw()




exp_10_diff_proportions_probs_percentiles_utilities %>% 
  arrange(desc(abs(median_right - median_left))) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity", fill = "dark gray") +
  geom_point(aes(y = p50, size = median_n_total_flies)) +
  geom_linerange(aes(ymin = p025, ymax = p975)) +
  theme_bw()




exp_10_diff_proportions_probs_percentiles_utilities %>% 
  arrange(desc(abs(median_right - median_left))) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_bar(aes(y = prop_right_more_than_left), stat = "identity", fill = "dark gray") +
  geom_point(aes(y = p50, size = median_n_total_flies)) +
  geom_linerange(aes(ymin = p025, ymax = p975)) +
  theme_bw()




exp_10_diff_proportions_probs_percentiles_utilities %>% 
  arrange(desc(abs(median_right - median_left))) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_point(aes(y = prop_right_more_than_left, 
                 size = median_n_total_flies),
             color = "dark gray") +
  geom_point(aes(y = p50), color = "black") +
  geom_linerange(aes(ymin = p025, ymax = p975), color = "black") +
  theme_bw()




exp_10_diff_proportions_probs_percentiles_utilities %>% 
  arrange(desc(abs(median_right - median_left))) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_point(aes(y = prop_right_more_than_left, 
                 size = median_n_total_flies),
             color = "dark gray") +
  geom_point(aes(y = p50), color = "black") +
  geom_linerange(aes(ymin = p025, ymax = p975), color = "black") +
  theme_bw()




exp_10_diff_proportions_probs_percentiles_utilities %>% 
  arrange(desc(median_rel_diff)) %>% 
  mutate(cs = 1:nrow(.)) %>% 
  ggplot(aes(x = cs)) +
  geom_point(aes(y = prop_right_more_than_left, 
                 size = median_n_total_flies),
             color = "dark gray") +
  geom_point(aes(y = p50), color = "black") +
  geom_linerange(aes(ymin = p025, ymax = p975), color = "black") +
  theme_bw()



View(exp_10_diff_proportions_probs_percentiles_utilities %>% mutate(abs_error = abs(p50 - prop_right_more_than_left)))
