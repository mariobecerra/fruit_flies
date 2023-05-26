library(tidyverse)
library(opdesmixr)
library(here)



i_opt_design_filename_flies_01 = paste0(here("european_fly/out/second_i_optimal_design.rds"))
i_opt_design_flies_01 = readRDS(i_opt_design_filename_flies_01)


which_min_i_opt = which.min(sapply(1:length(i_opt_design_flies_01), function(i) i_opt_design_flies_01[[i]]$opt_crit_value))

X_min = i_opt_design_flies_01[[which_min_i_opt]]$X

model_matrix_names = c("R", "O", "Y", "G", "B", "P", "R*O", "R*Y", "R*G", "R*B", "R*P", "R*UV", "O*Y", "O*G", "O*B", "O*P", "O*UV", "Y*G", "Y*B", "Y*P", "Y*UV", "G*B", "G*P", "G*UV", "B*P", "B*UV", "P*UV", "R*intensity", "O*intensity", "Y*intensity", "G*intensity", "B*intensity", "P*intensity", "UV*intensity", "intensity^2", "no_choice")

design_data_frame = mnl_design_array_to_dataframe(X_min) %>% 
  set_names(c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "cs"))



# Plot the real betas versus the draws
draws_beta_mean_sd = data.frame(
  variable = model_matrix_names,
  mean = apply(i_opt_design_flies_01[[which_min_i_opt]]$beta, 2, mean),
  sd = apply(i_opt_design_flies_01[[which_min_i_opt]]$beta, 2, sd)
) %>% 
  mutate(var_fct = fct_reorder(variable, sd))



draws_beta_mean_sd %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), size = 0.6, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1.2, position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")


draws_beta_mean_sd %>% 
  ggplot(aes(x = var_fct)) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), size = 0.6, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1.2, position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")




design_data_frame %>% 
  pivot_longer(R:intensity) %>% 
  ggplot() +
  geom_boxplot(aes(name, value))

design_data_frame %>% 
  pivot_longer(R:intensity) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~name, scales = "free")


SDs_prior = apply(i_opt_design_flies_01[[which_min_i_opt]]$beta, 2, sd)
names(SDs_prior) = model_matrix_names
sort(SDs_prior, decreasing = T)

IM_B = mnl_get_information_matrix(
  X_min, 
  beta = i_opt_design_flies_01[[which_min_i_opt]]$beta,
  order = 2,
  transform_beta = F,
  n_pv = 1,
  no_choice = T
  )
colnames(IM_B) = model_matrix_names
rownames(IM_B) = model_matrix_names

sort(diag(IM_B))





IM_avg = mnl_get_information_matrix(
  X_min, 
  beta = apply(i_opt_design_flies_01[[which_min_i_opt]]$beta, 2, mean),
  order = 2,
  transform_beta = F,
  n_pv = 1,
  no_choice = T
)
colnames(IM_avg) = model_matrix_names
rownames(IM_avg) = model_matrix_names

sort(diag(IM_avg))



write_csv(design_data_frame, here("european_fly/out/second_i_optimal_design.csv"))



model_matrix = matrix(0.0, ncol = mnl_get_number_of_parameters(dim(X_min)[1]-1, 1, 2)-1, nrow = dim(X_min)[2]*dim(X_min)[3])
colnames(model_matrix) = model_matrix_names[1:(length(model_matrix_names)-1)]
for(s in 1:dim(X_min)[3]){
  init = dim(X_min)[2]*(s - 1) + 1
  end = dim(X_min)[2]*s
  model_matrix[init:end, ] = mnl_get_Xs(X_min, s, 2, 1, T)
}

model_matrix %>% 
  as.data.frame() %>% 
  mutate(ix = 1:nrow(.)) %>% 
  pivot_longer(-ix) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~name, scales = "free")





