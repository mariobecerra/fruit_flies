# Based on what Ya sent.
# To transform from the mixture space to the format that Frank uses.

library(readxl)
library(tidyverse)
library(opdesmixr)
library(here)

source(here("japanese_fly/code/transform_design_function.R"))

# IMPORT DESIGN --------------------------------------------------------------
intensity_df <- read_excel(here("japanese_fly/data/intensities_thorlabs.xlsx")) %>% 
  as.data.frame()

design_object = readRDS(here("japanese_fly/out/japanese_flies_third_i_optimal_design.rds"))

i_opts = sapply(seq_along(design_object), function(i) design_object[[i]]$opt_crit_value)

design_df = mnl_design_array_to_dataframe(design_object[[which.min(i_opts)]]$X) %>% 
  set_names(c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "cs")) %>% 
  mutate(cs = as.numeric(cs)) %>% 
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


