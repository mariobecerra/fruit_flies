library(tidyverse)
library(here)

design_1_mapping = read_csv(here("data/01 DMela_count_per_image.csv")) %>% 
  select(folder, 
         R_left, 
         O_left, 
         Y_left, 
         G_left, 
         B_left, 
         P_left, 
         UV_left, 
         intensity_left, 
         R_right, 
         O_right, 
         Y_right, 
         G_right, 
         B_right, 
         P_right, 
         UV_right, 
         intensity_right
  ) %>% 
  distinct()





design_same_colors = read_csv(here("european_fly/out/design_same_colors_UV_w_intensity.csv")) %>% 
  group_by(choice_set) %>% 
  mutate(alternative = 1:n()) %>% 
  ungroup()


design_same_colors_transformed_01 = read_tsv(here("data/design_same_colors_UV_transformed.txt"))


key_design_same_colors_df_01 = design_same_colors_transformed_01 %>% 
  select(5:ncol(.)) %>% 
  mutate_at(vars(starts_with("Var")), ~as.character(round(.))) %>% 
  mutate_at(vars(!starts_with("Var")), ~substr(as.character(.), 1, 1)) %>% 
  distinct()

keys_design_same_colors_01 = apply(key_design_same_colors_df_01, 1, function(x) paste0(x, collapse = ""))


names_colors_01 = design_same_colors_transformed_01 %>% 
  select(5:ncol(.)) %>% 
  select(!starts_with("Var")) %>% 
  names()

# design_same_colors_transformed_2 = design_same_colors_transformed %>% 
#   select(starts_with("Var")) %>% 
#   set_names(c(paste0(names_colors[1:(length(names_colors)/2)], "_left"), 
#               paste0(names_colors[1:(length(names_colors)/2)], "_right"))) %>% 
#   distinct() %>% 
#   mutate(key = keys_design_same_colors) %>% 
#   distinct()



design_same_colors_mapping_01 = design_same_colors %>% 
  filter(alternative == 1) %>% 
  set_names(paste0(names(design_same_colors), "_left")) %>% 
  bind_cols(
    design_same_colors %>% 
      filter(alternative == 2) %>% 
      set_names(paste0(names(design_same_colors), "_right"))
  ) %>% 
  select(-starts_with("choice_set"), -starts_with("alternative")) %>% 
  mutate(folder = keys_design_same_colors_01)






design_same_colors_transformed_02 = read_tsv(here("data/design_same_colors_UV_transformed_2.txt"))


key_design_same_colors_df_02 = design_same_colors_transformed_02 %>% 
  select(5:ncol(.)) %>% 
  mutate_at(vars(starts_with("Var")), ~as.character(round(.))) %>% 
  mutate_at(vars(!starts_with("Var")), ~substr(as.character(.), 1, 1)) %>% 
  distinct()

keys_design_same_colors_02 = apply(key_design_same_colors_df_02, 1, function(x) paste0(x, collapse = ""))


names_colors_02 = design_same_colors_transformed_02 %>% 
  select(5:ncol(.)) %>% 
  select(!starts_with("Var")) %>% 
  names()


design_same_colors_mapping_02 = design_same_colors %>% 
  filter(alternative == 1) %>% 
  set_names(paste0(names(design_same_colors), "_left")) %>% 
  bind_cols(
    design_same_colors %>% 
      filter(alternative == 2) %>% 
      set_names(paste0(names(design_same_colors), "_right"))
  ) %>% 
  select(-starts_with("choice_set"), -starts_with("alternative")) %>% 
  mutate(folder = keys_design_same_colors_02)












design_1_swapped_mapping = read_csv(here("data/mapping_idesign_swapped.csv")) %>% 
  select(-chid)

design_1_swapped_mapping2 = read_csv(here("data/mapping_idesign_swapped2.csv"))















second_i_opt_design_european = read_csv(here("european_fly/out/second_i_optimal_design.csv")) %>% 
  rename(choice_set = cs) %>% 
  group_by(choice_set) %>% 
  mutate(alternative = 1:n()) %>% 
  ungroup()



second_i_opt_design_european_transformed_01 = read_tsv(here("data/second_i_optimal_design_transformed.txt"))


key_second_i_opt_design_european_df_01 = second_i_opt_design_european_transformed_01 %>% 
  filter(Image != 300) %>% 
  select(5:ncol(.)) %>% 
  mutate_at(vars(starts_with("Var")), ~as.character(round(.))) %>% 
  mutate_at(vars(!starts_with("Var")), ~substr(as.character(.), 1, 1))

keys_second_i_opt_design_european_01 = apply(key_second_i_opt_design_european_df_01, 1, function(x) paste0(x, collapse = ""))



second_i_opt_design_european_mapping_01 = second_i_opt_design_european %>% 
  filter(alternative == 1) %>% 
  set_names(paste0(names(second_i_opt_design_european), "_left")) %>% 
  bind_cols(
    second_i_opt_design_european %>% 
      filter(alternative == 2) %>% 
      set_names(paste0(names(second_i_opt_design_european), "_right"))
  ) %>% 
  select(-starts_with("choice_set"), -starts_with("alternative")) %>% 
  mutate(folder = keys_second_i_opt_design_european_01)










# Not really swapped, Ya made a mistake and reversed the order of the choice sets instead of the alternatives within each choice set.
# second_i_opt_design_european_swapped = second_i_opt_design_european %>% 
#   mutate(alternative = ifelse(alternative == 1, 2, 1)) %>% 
#   arrange(choice_set, alternative)









design_mapping = design_1_mapping %>% 
  bind_rows(design_same_colors_mapping_01) %>% 
  bind_rows(design_same_colors_mapping_02) %>% 
  bind_rows(design_1_swapped_mapping) %>% 
  bind_rows(design_1_swapped_mapping2) %>% 
  bind_rows(second_i_opt_design_european_mapping_01)



write_csv(design_mapping, here("european_fly/out/design_mapping_european.csv"))





