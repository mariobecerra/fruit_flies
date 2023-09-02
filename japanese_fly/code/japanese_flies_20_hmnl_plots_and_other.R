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


########################################
##### Parameter plots
########################################

names_betas_level_0 = c("R", "O", "Y", "G", "B", "P", "R*O", "R*Y", "R*G", "R*B", "R*P", "R*UV", "O*Y", "O*G", "O*B", "O*P", "O*UV", "Y*G", "Y*B", "Y*P", "Y*UV", "G*B", "G*P", "G*UV", "B*P", "B*UV", "P*UV", "R*intensity", "O*intensity", "Y*intensity", "G*intensity", "B*intensity", "P*intensity", "UV*intensity", "intensity^2", "no_choice")

betas_level_0_summary = summary(model_stan, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = names_betas_level_0,
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))



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




hyper_betas_plot = betas_level_0_summary %>%
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(y = mean), color = "red", size = 1.5) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "red", size = 0.5) +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "red", size = 1.3) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle(paste0("Japanese fly hyper betas, experiments 1 to ", n_experiments))


ggsave(plot = hyper_betas_plot,
       filename = here("japanese_fly/out/plots/parameter_hyperbetas.png"),
       width = 18, height = 12, units = "cm", dpi = 300)


ggsave(plot = hyper_betas_plot + ggtitle("") + ylab("") + xlab(""),
       filename = here("japanese_fly/out/plots/parameter_hyperbetas_2.png"),
       width = 15, height = 12, units = "cm", dpi = 150)

knitr::plot_crop(here("japanese_fly/out/plots/parameter_hyperbetas_2.png"))


plot_all_betas_gray_red_01 = betas_level_0_summary %>% 
  mutate(exp = "All experiments", line_thickness_1 = 2, line_thickness_2 = 0.8) %>% 
  bind_rows(betas_level_1_summary %>% 
              mutate(line_thickness_1 = 0.8, line_thickness_2 = 0.5)) %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.9), size = 0.6) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd, size = line_thickness_2), 
                 position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd, size = line_thickness_1),
                 position = position_dodge(width = 0.9), show.legend = FALSE) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") +
  ggtitle(paste0("Japanese fly betas of experiments 1 to ", n_experiments)) +
  scale_size_continuous(range = c(0.7, 1.2)) +
  scale_color_manual(values = c("red", rep("dark grey", n_experiments)))


ggsave(plot = plot_all_betas_gray_red_01,
       filename = here("japanese_fly/out/plots/parameter_all_betas_gray_red_01.png"),
       width = 30, height = 15, units = "cm", dpi = 300)


plot_all_betas_gray_red_02 = betas_level_0_summary %>% 
  mutate(exp = "All experiments", line_thickness_2 = 0.8) %>% 
  bind_rows(betas_level_1_summary %>% 
              mutate(line_thickness_2 = 0.5)) %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.4) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.9), size = 0.6) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd, size = line_thickness_2), 
                 position = position_dodge(width = 0.9), show.legend = FALSE) +
  # geom_linerange(aes(ymin = mean - sd, ymax = mean + sd, size = line_thickness_1),
  #                position = position_dodge(width = 0.7), show.legend = FALSE) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.4, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") +
  ggtitle(paste0("Japanese fly betas of experiments 1 to ", n_experiments)) +
  scale_size_continuous(range = c(0.7, 1.1)) +
  scale_color_manual(values = c("red", rep("dark grey", n_experiments)))


ggsave(plot = plot_all_betas_gray_red_02,
       filename = here("japanese_fly/out/plots/parameter_all_betas_gray_red_02.png"),
       width = 30, height = 15, units = "cm", dpi = 300)


ggsave(plot = plot_all_betas_gray_red_02 + ggtitle("") + ylab("") + xlab(""),
       filename = here("japanese_fly/out/plots/parameter_all_betas_gray_red_03.png"),
       width = 22, height = 12, units = "cm", dpi = 150)

knitr::plot_crop(here("japanese_fly/out/plots/parameter_all_betas_gray_red_03.png"))




Sigma_level_0_posterior_median = apply(Sigma_level_0_draws_hmnl$Sigma_level_0, 2:3, median)
Sigma_level_0_posterior_p025 = apply(Sigma_level_0_draws_hmnl$Sigma_level_0, 2:3, quantile, 0.025)
Sigma_level_0_posterior_p10 = apply(Sigma_level_0_draws_hmnl$Sigma_level_0, 2:3, quantile, 0.1)
Sigma_level_0_posterior_p90 = apply(Sigma_level_0_draws_hmnl$Sigma_level_0, 2:3, quantile, 0.9)
Sigma_level_0_posterior_p975 = apply(Sigma_level_0_draws_hmnl$Sigma_level_0, 2:3, quantile, 0.975)

Sigma_level_0_posterior_diag_df = data.frame(
  variable = names_betas_level_0,
  p025_sd = sqrt(diag(Sigma_level_0_posterior_p025)),
  p10_sd = sqrt(diag(Sigma_level_0_posterior_p10)),
  median_sd = sqrt(diag(Sigma_level_0_posterior_median)),
  p90_sd = sqrt(diag(Sigma_level_0_posterior_p90)),
  p975_sd = sqrt(diag(Sigma_level_0_posterior_p975))
  )

Sigma_level_0_posterior_diag_df

write_csv(Sigma_level_0_posterior_diag_df, 
          here("japanese_fly/out/plots/Sigma_level_0_posterior_diag_df.csv"))



Sigma_level_0_posterior_diag_df_plot = Sigma_level_0_posterior_diag_df %>% 
  mutate(ix = 1:nrow(.)) %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable)) +
  xlab("") +
  ylab("") +
  geom_point(aes(y = median_sd), size = 1.5) +
  geom_linerange(aes(ymin = p025_sd, ymax = p975_sd), size = 0.5) +
  geom_linerange(aes(ymin = p10_sd, ymax = p90_sd), size = 1.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ggtitle("Standard deviation of beta")


ggsave(plot = Sigma_level_0_posterior_diag_df_plot,
       filename = here("japanese_fly/out/plots/Sigma_level_0_posterior_diag.png"),
       width = 15, height = 9, units = "cm", dpi = 150)


ggsave(plot = Sigma_level_0_posterior_diag_df_plot + ggtitle(""),
       filename = here("japanese_fly/out/plots/Sigma_level_0_posterior_diag_2.png"),
       width = 15, height = 9, units = "cm", dpi = 150)

knitr::plot_crop(here("japanese_fly/out/plots/Sigma_level_0_posterior_diag_2.png"))







########################################
##### Data plots
########################################

design_5_with_mixture = read.csv(here("japanese_fly/out/5th_i_optimal_design_japanese_mapping.csv"))



# The first 22 rows here are the 22 choice sets that I created and forcefully added in the experiment
# The last 6 rows have one fixed alternative of the 6 handpicked colors I wanted to check out
hand_picked_choice_sets = design_5_with_mixture[1:22,]



# hand_picked_choice_sets_tibble = readRDS(here("japanese_fly/out/handpicked_choice_sets_5th_i_opt_design.rds"))
# 
# # filtered by hand the choice sets with very very similar alternatives
# # This is because the coordinate exchange was having trouble computing the I-optimality
# hand_picked_choice_sets_tibble_filtered = hand_picked_choice_sets_tibble %>% 
#   filter(!(choice_set %in% c(29, 30, 3, 4, 1, 2, 11, 12))) %>% 
#   select(-alt)


hand_picked_colors = readRDS(here("japanese_fly/out/handpicked_colors_5th_i_opt_design.rds"))



choice_sets_to_analyze = design_5_with_mixture[1:22, "folder"]
# choice_sets_to_analyze = design_5_with_mixture[1:28, "folder"]





filtered_exp_9_cs_of_interest = counts_japanese %>% 
  filter(experiment == 9, folder %in% choice_sets_to_analyze) %>% 
  select(R:no_choice)

filtered_exp_9_cs_of_interest_to_predict_list <- create_model_matrix_second_order_scheffe(filtered_exp_9_cs_of_interest)


predicted_point_utilities = filtered_exp_9_cs_of_interest_to_predict_list$X %*% betas_level_0_summary$mean


filtered_exp_9_cs_of_interest_utilities = counts_japanese %>% 
  filter(experiment == 9, folder %in% choice_sets_to_analyze) %>% 
  mutate(pred_ut = as.numeric(predicted_point_utilities))



filtered_exp_9_cs_of_interest_utilities

filtered_exp_9_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") 




choice_sets_of_interest_catalog = filtered_exp_9_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  select(folder, choice_set, alternative, pred_ut, R:intensity) %>% 
  distinct() %>% 
  mutate(color = paste(R, O, Y, G, B, P, UV, intensity)) %>% 
  mutate(color_ix = as.integer(as.factor(color)))


exp_9_handpicked_colors_indices_and_utilities = choice_sets_of_interest_catalog %>% 
  mutate(exp_pred_ut = exp(pred_ut)) %>% 
  select(color_ix, pred_ut, exp_pred_ut, R:intensity) %>% 
  distinct() %>% 
  arrange(color_ix)

exp_9_handpicked_colors_indices_and_utilities

write_csv(exp_9_handpicked_colors_indices_and_utilities, 
          here("japanese_fly/out/plots/exp_9_handpicked_colors_indices_and_utilities.csv"))



width_counts_plots = 32
height_counts_plots = 20


counts_exp9_handpicked_colors_01 = filtered_exp_9_cs_of_interest_utilities %>% 
  filter(alternative != "3_none") %>% 
  ggplot() +
  geom_jitter(aes(x = alternative, y = count), size = 0.6) +
  facet_wrap(~folder) +
  theme_bw()


ggsave(
  plot = counts_exp9_handpicked_colors_01,
  filename = here("japanese_fly/out/plots/counts_exp9_handpicked_colors_01.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)




utilities_counts_tibble = filtered_exp_9_cs_of_interest_utilities %>% 
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



counts_exp9_handpicked_colors_02 = utilities_counts_tibble %>% 
  ggplot(aes(x = alternative, y = median_counts)) +
  geom_point() +
  geom_linerange(aes(ymin = p10_counts, ymax = p90_counts)) +
  facet_wrap(~group) +
  theme_bw() +
  ggtitle("Percentiles 10 and 90 of counts in each of the choice sets of interest")


ggsave(
  plot = counts_exp9_handpicked_colors_02,
  filename = here("japanese_fly/out/plots/counts_exp9_handpicked_colors_02.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)



counts_exp9_handpicked_colors_03 = filtered_exp_9_cs_of_interest_utilities %>% 
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
  plot = counts_exp9_handpicked_colors_03,
  filename = here("japanese_fly/out/plots/counts_exp9_handpicked_colors_03.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)



counts_exp9_handpicked_colors_04 = filtered_exp_9_cs_of_interest_utilities %>% 
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
  theme_bw()

ggsave(
  plot = counts_exp9_handpicked_colors_04,
  filename = here("japanese_fly/out/plots/counts_exp9_handpicked_colors_04.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)




counts_exp9_handpicked_colors_05 = filtered_exp_9_cs_of_interest_utilities %>% 
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
  theme_bw()


ggsave(
  plot = counts_exp9_handpicked_colors_05,
  filename = here("japanese_fly/out/plots/counts_exp9_handpicked_colors_05.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)





counts_exp9_handpicked_colors_06 = filtered_exp_9_cs_of_interest_utilities %>% 
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
  theme_bw()


ggsave(
  plot = counts_exp9_handpicked_colors_06,
  filename = here("japanese_fly/out/plots/counts_exp9_handpicked_colors_06.png"),
  width = width_counts_plots, height = height_counts_plots, units = "cm", dpi = 300
)

## Comparisons missing and that should perhaps be there
# 25, 35, 34










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

