library(rstan)
library(MASS)
library(tidyverse)
library(ggtern)
library(here)

source(here("japanese_fly/code/japanese_flies_hmnl_1_level_utils.R"))


model_stan_1_to_4 = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_4_stan_object_30images.rds"))
model_stan_1_to_7 = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_7_stan_object_30images.rds"))
model_stan_1_to_8 = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_8_stan_object_30images.rds"))
model_stan_1_to_9 = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_9_stan_object_30images.rds"))




counts_japanese = read_csv(here("japanese_fly/out/counts_japanese_choice_format.csv"))

hand_picked_choice_sets_4th_design = readRDS(here("japanese_fly/out/handpicked_choice_sets_4th_i_opt_design.rds")) %>% 
  select(-alt) %>% 
  rename(cs = choice_set)



hand_picked_choice_sets_5th_design = readRDS(here("japanese_fly/out/handpicked_choice_sets_5th_i_opt_design.rds")) %>% 
  select(-alt) %>% 
  rename(cs = choice_set)


hand_picked_choice_sets_folders_4th_design = c(
  "T2F0F0T12F0F0T18F0T7F0F0F0F0F0F0F0",
  "T2F0F0T12F0F0T18F0T1F0F0T11F0F0T17F0",
  "T13F0F0F0F0F0F0F0T1F0F0T11F0F0T18F0",
  "T2F0F0T12F0F0T17F0T1F0F0T11F0F0T18F0",
  "T2F0F0T12F0F0T18F0T1F0F0T12F0F0T18F0",
  "T2F0F0T13F0F0T19F0T1F0F0T11F0F0T18F0",
  "T13F0F0F0F0F0F0F0T1F0F0T12F0F0T18F0",
  "T13F0F0F0F0F0F0F0T1F0F0T11F0F0T17F0",
  "T2F0F0T13F0F0T19F0T1F0F0T11F0F0T17F0",
  "T2F0F0T12F0F0T17F0T7F0F0F0F0F0F0F0",
  "T2F0F0T13F0F0T19F0T7F0F0F0F0F0F0F0",
  "T2F0F0T12F0F0T17F0T1F0F0T12F0F0T18F0")





library(readxl)

intensity_df <- read_excel(here("japanese_fly/data/intensities_thorlabs.xlsx")) %>% 
  as.data.frame()

source(here("japanese_fly/code/transform_design_function.R"))

hand_picked_choice_sets_folders_5th_design =  transform_design_flies(
  design_df = hand_picked_choice_sets_5th_design, 
  intensity_df = intensity_df
)$design_with_mixture$folder





counts_japanese %>% 
  filter(experiment == 8,
         folder %in% hand_picked_choice_sets_folders_4th_design)


counts_japanese %>% 
  filter(experiment == 9,
         folder %in% hand_picked_choice_sets_folders_5th_design)



counts_handpicked_cs_exp8 = counts_japanese %>% 
  filter(experiment == 8,
         folder %in% hand_picked_choice_sets_folders_4th_design) %>% 
  # mutate(color = paste(R, O, Y, G, B, P, UV, intensity, no_choice, sep = ", ")) %>% 
  mutate(color = paste(round(R, 3), round(O, 3), round(Y, 3), round(G, 3), round(B, 3), round(P, 3), round(UV, 3), round(intensity, 3), no_choice, sep = ", ")) %>% 
  mutate(color_int = as.integer(as.factor(color)))


counts_handpicked_cs_exp8 %>% 
  filter(no_choice != 1) %>%
  group_by(color, color_int) %>% 
  summarize(
    min_flies = min(count),
    median_flies = median(count),
    mean_flies = mean(count),
    max_flies = max(count),
  )


counts_handpicked_cs_exp8 %>% 
  filter(no_choice != 1) %>%
  ggplot() +
  geom_histogram(aes(count)) +
  facet_wrap(~color, scales = "free_y")
  





counts_handpicked_cs_exp9 = counts_japanese %>% 
  filter(experiment == 9,
         folder %in% hand_picked_choice_sets_folders_5th_design) %>% 
  # mutate(color = paste(R, O, Y, G, B, P, UV, intensity, no_choice, sep = ", ")) %>% 
  mutate(color = paste(round(R, 3), round(O, 3), round(Y, 3), round(G, 3), round(B, 3), round(P, 3), round(UV, 3), round(intensity, 3), no_choice, sep = ", ")) %>% 
  mutate(color_int = as.integer(as.factor(color)))


counts_handpicked_cs_exp9 %>% 
  filter(no_choice != 1) %>%
  group_by(color, color_int) %>% 
  summarize(
    min_flies = min(count),
    median_flies = median(count),
    mean_flies = mean(count),
    max_flies = max(count),
  )


counts_handpicked_cs_exp9 %>% 
  filter(no_choice != 1) %>%
  ggplot() +
  geom_histogram(aes(count)) +
  facet_wrap(~color, scales = "free_y")





counts_japanese %>% 
  filter(experiment == 8) %>% 
  # mutate(color = paste(R, O, Y, G, B, P, UV, intensity, no_choice, sep = ", ")) %>% 
  mutate(color = paste(round(R, 3), round(O, 3), round(Y, 3), round(G, 3), round(B, 3), round(P, 3), round(UV, 3), round(intensity, 3), no_choice, sep = ", ")) %>% 
  mutate(color_int = as.integer(as.factor(color))) %>% 
  filter(no_choice != 1) %>%
  group_by(color, color_int) %>% 
  summarize(
    min_flies = min(count),
    median_flies = median(count),
    mean_flies = mean(count),
    max_flies = max(count),
  ) %>% 
  arrange(desc(mean_flies))





counts_japanese %>% 
  filter(experiment == 9) %>% 
  # mutate(color = paste(R, O, Y, G, B, P, UV, intensity, no_choice, sep = ", ")) %>% 
  mutate(color = paste(round(R, 3), round(O, 3), round(Y, 3), round(G, 3), round(B, 3), round(P, 3), round(UV, 3), round(intensity, 3), no_choice, sep = ", ")) %>% 
  mutate(color_int = as.integer(as.factor(color))) %>% 
  filter(no_choice != 1) %>%
  group_by(color, color_int) %>% 
  summarize(
    min_flies = min(count),
    median_flies = median(count),
    mean_flies = mean(count),
    max_flies = max(count),
  ) %>% 
  arrange(desc(mean_flies))










names_betas_level_0 = c("R", "O", "Y", "G", "B", "P", "R*O", "R*Y", "R*G", "R*B", "R*P", "R*UV", "O*Y", "O*G", "O*B", "O*P", "O*UV", "Y*G", "Y*B", "Y*P", "Y*UV", "G*B", "G*P", "G*UV", "B*P", "B*UV", "P*UV", "R*intensity", "O*intensity", "Y*intensity", "G*intensity", "B*intensity", "P*intensity", "UV*intensity", "intensity^2", "no_choice")


betas_level_0_summary_1_to_4 = summary(model_stan_1_to_4, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = names_betas_level_0,
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))


betas_level_0_summary_1_to_7 = summary(model_stan_1_to_7, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = names_betas_level_0,
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))


betas_level_0_summary_1_to_8 = summary(model_stan_1_to_8, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = names_betas_level_0,
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))


betas_level_0_summary_1_to_9 = summary(model_stan_1_to_9, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = names_betas_level_0,
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))




betas_level_0_summary_1_to_4 %>% 
  mutate(model = "Exp 1 to 4") %>% 
  bind_rows(
    betas_level_0_summary_1_to_7 %>% 
      mutate(model = "Exp 1 to 7")
  ) %>% 
  bind_rows(
    betas_level_0_summary_1_to_8 %>% 
      mutate(model = "Exp 1 to 8")
  ) %>% 
  bind_rows(
    betas_level_0_summary_1_to_9 %>% 
      mutate(model = "Exp 1 to 9")
  ) %>% 
  ggplot(aes(x = variable, color = model)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.8)) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = 0.8)) +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), size = 1, position = position_dodge(width = 0.8)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Japanese fly hyper betas")





