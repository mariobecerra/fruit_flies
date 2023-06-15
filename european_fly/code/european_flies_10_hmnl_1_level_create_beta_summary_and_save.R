library(here)
library(tidyverse)
library(rstan)


model_stan = readRDS(here("european_fly/out/09_hmnl_1_level_03_stan_object_5000_iter_50_images.rds"))

names_betas_level_0 = c("R", "O", "Y", "G", "B", "P", "R*O", "R*Y", "R*G", "R*B", "R*P", "R*UV", "O*Y", "O*G", "O*B", "O*P", "O*UV", "Y*G", "Y*B", "Y*P", "Y*UV", "G*B", "G*P", "G*UV", "B*P", "B*UV", "P*UV", "R*intensity", "O*intensity", "Y*intensity", "G*intensity", "B*intensity", "P*intensity", "UV*intensity", "intensity^2", "no_choice")


betas_level_0_summary = summary(model_stan, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = names_betas_level_0,
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = T))


saveRDS(betas_level_0_summary, here("european_fly/out/09_hmnl_1_level_03_betas_level_0_summary.rds"))

betas_level_0_summary %>% 
  select(variable, mean, sd) %>% 
  write_csv(., here("european_fly/out/09_hmnl_1_level_03_betas_level_0_mean_sd.csv"))


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




