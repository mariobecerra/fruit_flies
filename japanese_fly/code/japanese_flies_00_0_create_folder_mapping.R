library(tidyverse)
library(here)

design_1_mapping = read_csv(here("japanese_fly/data/design_1_mapping.csv"))

design_2_mapping = read_csv(here("japanese_fly/out/2nd_i_optimal_design_japanese_mapping.csv"))

design_3_mapping = read_csv(here("japanese_fly/out/3rd_i_optimal_design_japanese_mapping.csv"))

# design_4_mapping = read_csv(here("japanese_fly/out/4th_i_optimal_design_japanese_mapping.csv"))


design_mapping = design_1_mapping %>%
  bind_rows(design_2_mapping) %>% 
  bind_rows(design_3_mapping)

write_csv(design_mapping, "japanese_fly/out/design_mapping_japanese.csv")





