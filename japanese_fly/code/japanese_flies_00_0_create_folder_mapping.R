library(tidyverse)
library(here)

design_1_mapping = read_csv(here("japanese_fly/data/design_1_mapping.csv"))

design_2_mapping = read_csv(here("japanese_fly/out/2nd_i_optimal_design_japanese_mapping.csv"))

design_3_mapping = read_csv(here("japanese_fly/out/3rd_i_optimal_design_japanese_mapping.csv"))

# In case of design 4, there were two points in the design that were very similar, so the folder name is the same. We'll remove one of them because for practical purposes we can use any of the two points and the parameters wouldn't change
design_4_mapping = read_csv(here("japanese_fly/out/4th_i_optimal_design_japanese_mapping.csv")) %>% 
  filter(!duplicated(folder))


design_mapping = design_1_mapping %>%
  bind_rows(design_2_mapping) %>% 
  bind_rows(design_3_mapping) %>% 
  bind_rows(design_4_mapping)

write_csv(design_mapping, "japanese_fly/out/design_mapping_japanese.csv")





