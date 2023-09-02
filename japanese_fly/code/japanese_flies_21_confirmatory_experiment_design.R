library(tidyverse)
library(readxl)
library(opdesmixr)
library(here)

var_names = c("R", "O", "Y", "G", "B", "P", "R*O", "R*Y", "R*G", "R*B", "R*P", "R*UV", "O*Y", "O*G", "O*B", "O*P", "O*UV", "Y*G", "Y*B", "Y*P", "Y*UV", "G*B", "G*P", "G*UV", "B*P", "B*UV", "P*UV", "R*intensity", "O*intensity", "Y*intensity", "G*intensity", "B*intensity", "P*intensity", "UV*intensity", "intensity^2", "no_choice")


q_flies_01 = 7
J_flies_01 = 2
S_flies_01 = 80
order_flies_01 = 2
n_pv_flies_01 = 1
no_choice_flies_01 = T

# read hand picked choice sets
hand_picked_choice_sets_tibble = readRDS(here("japanese_fly/out/handpicked_choice_sets_5th_i_opt_design.rds"))




# Based on what Ya sent.
# To transform from the mixture space to the format that Frank uses.


source(here("japanese_fly/code/transform_design_function.R"))

intensity_df <- read_excel(here("japanese_fly/data/intensities_thorlabs.xlsx")) %>% 
  as.data.frame()

# Read last design
design_object = readRDS(here("japanese_fly/out/japanese_flies_5th_i_optimal_design.rds"))

i_opts = sapply(seq_along(design_object), function(i) design_object[[i]]$opt_crit_value)

design_5th_df = mnl_design_array_to_dataframe(design_object[[which.min(i_opts)]]$X) %>% 
  set_names(c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "cs")) %>% 
  mutate(cs = as.numeric(cs)) %>% 
  as.data.frame()


# Keep only the 23rd choice set onwards because the first 22 are the ones that were
# manually set and manually filtered because of the computation problem in the information matrix
design_5th_df_filtered = design_5th_df %>% 
  filter(cs >= 23) %>% 
  mutate(choice_set = paste0("d5_cs", as.integer(as.factor(cs)))) %>% 
  select(-cs)


# The first 30 choice sets are the hand picked ones
# The next 6 are the ones with the handpickeed colors and another collor chosen by the coord exchange
design_df_ordered_1 = hand_picked_choice_sets_tibble %>% 
  mutate(choice_set = paste0("handpicked_cs", as.integer(as.factor(choice_set)))) %>% 
  select(-alt) %>% 
  bind_rows(design_5th_df_filtered) %>% 
  mutate(cs = rep(1:(nrow(.)/J_flies_01), each = J_flies_01),
         alt = rep(1:J_flies_01, times = nrow(.)/2))

set.seed(2023)
cs_order_ix = tibble(
  cs_original = sample(1:max(design_df_ordered_1$cs), max(design_df_ordered_1$cs))
) %>% 
  mutate(cs_new = 1:nrow(.))


design_df_ordered_2 = design_df_ordered_1 %>% 
  inner_join(
    cs_order_ix, by = c("cs" = "cs_original")
  ) 

write.csv(design_df_ordered_2,
          here("japanese_fly/out/confirmatory_design_japanese_ordered.csv"), 
          row.names = F)

design_df = design_df_ordered_2 %>% 
  select(-cs, -choice_set) %>% 
  rename(cs = cs_new) %>% 
  arrange(cs, alt) %>% 
  select(-alt)

transformed_design_list <- transform_design_flies(
  design_df = design_df, 
  intensity_df = intensity_df
)




write.csv(transformed_design_list$design_software, 
          here("japanese_fly/out/confirmatory_design_japanese_frank.csv"), 
          row.names = F)

write.csv(transformed_design_list$design_with_mixture, 
          here("japanese_fly/out/confirmatory_design_japanese_mapping.csv"), 
          row.names = F)



