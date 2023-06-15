# To transform from the mixture space to the format that Frank uses.
# Also, add a few hand-chosen choice sets to confirm the favorite color for the flies.

library(readxl)
library(tidyverse)
library(opdesmixr)
library(here)

source(here("japanese_fly/code/transform_design_function.R"))

# read hand picked choice sets
hand_picked_choice_sets = readRDS(here("japanese_fly/out/handpicked_choice_sets_4th_i_opt_design.rds")) %>% 
  select(-alt) %>% 
  rename(cs = choice_set)

# IMPORT DESIGN --------------------------------------------------------------
intensity_df <- read_excel(here("japanese_fly/data/intensities_thorlabs.xlsx")) %>% 
  as.data.frame()

design_object = readRDS(here("japanese_fly/out/japanese_flies_4th_i_optimal_design.rds"))

i_opts = sapply(seq_along(design_object), function(i) design_object[[i]]$opt_crit_value)

design_df_0 = mnl_design_array_to_dataframe(design_object[[which.min(i_opts)]]$X) %>% 
  set_names(c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "cs")) %>% 
  mutate(cs = as.numeric(cs) + max(hand_picked_choice_sets$cs)) %>% 
  bind_rows(hand_picked_choice_sets)

set.seed(2023)
cs_order = sample(unique(design_df_0$cs))

design_df = tibble(cs = cs_order) %>% 
  left_join(design_df_0) %>% 
  select(-cs) %>% 
  mutate(cs = rep(1:length(cs_order), each = 2)) %>% 
  as.data.frame()
  



transformed_design_list <- transform_design_flies(
  design_df = design_df, 
  intensity_df = intensity_df
)




write.csv(transformed_design_list$design_software, here("japanese_fly/out/4th_i_optimal_design_japanese_frank.csv"), row.names = F)
write.csv(transformed_design_list$design_with_mixture, here("japanese_fly/out/4th_i_optimal_design_japanese_mapping.csv"), row.names = F)



design_df %>%
  pivot_longer(1:8) %>%
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~name, scales = "free_y")


