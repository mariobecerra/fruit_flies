# Based on what Ya sent.
# To transform from the mixture space to the format that Frank uses.

library(readxl)
library(tidyverse)
library(opdesmixr)
library(here)

source(here("code/transform_design_function.R"))

# IMPORT DESIGN --------------------------------------------------------------
intensity_df <- read_excel(here("data/intensity.xlsx")) %>% 
  as.data.frame()

# design_object = readRDS(here("out/second_i_optimal_design_japanese_old.rds"))
design_object = readRDS(here("out/second_i_optimal_design_japanese.rds"))

i_opts = sapply(seq_along(design_object), function(i) design_object[[i]]$opt_crit_value)

design_df = mnl_design_array_to_dataframe(design_object[[which.min(i_opts)]]$X) %>% 
  set_names(c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "cs")) %>% 
  mutate(cs = as.numeric(cs)) %>% 
  as.data.frame()



transformed_design_list <- transform_design_flies(
  design_df = design_df, 
  intensity_df = intensity_df
  )



write_csv(transformed_design_list$design_software, here("out/2nd_i_optimal_design_japanese_frank.csv"))
write_csv(transformed_design_list$design_with_mixture, here("out/2nd_i_optimal_design_japanese_mapping.csv"))



