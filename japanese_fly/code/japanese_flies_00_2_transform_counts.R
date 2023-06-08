library(lubridate)
library(tidyverse)
library(here)


# Update with latest counts!!!
count_data = read_csv(here('japanese_fly/out/counts_japanese_2023-06-08.csv'))


design_mapping = read_csv(here("japanese_fly/out/design_mapping_japanese.csv"))

# Assume 80 flies per experiment
n_experiments = max(count_data$experiment)
n_flies_per_experiment = tibble(
  experiment = 1:n_experiments,
  n_total_flies = rep(80, n_experiments)
)

merged_counts = count_data %>% 
  mutate(date_time = ymd_hms(paste(date, time_of_day))) %>% 
  select(experiment, date_time, choice_set, folder, count_left, count_right) %>% 
  left_join(design_mapping, by = "folder") %>% 
  left_join(n_flies_per_experiment, by = "experiment") %>% 
  arrange(date_time) %>% 
  mutate(count_no_choice = n_total_flies - count_left - count_right)

vars_to_select = c("experiment", "date_time", "folder", "choice_set")

names_mapping = count_data %>% 
  select(ends_with("left")) %>% 
  names() %>% 
  as.data.frame() %>% 
  set_names("original_name_left") %>% 
  mutate(
    new_name = stringi::stri_replace_first(str = original_name_left, fixed = "_left", "")
  ) %>% 
  mutate(new_name = ifelse(new_name == "UV1", "UV", ifelse(new_name == "UV2", "intensity", new_name)))



merged_left = merged_counts %>% 
  select(all_of(vars_to_select), ends_with("left")) %>% 
  set_names(c(vars_to_select, names_mapping$new_name)) %>% 
  mutate(
    alternative = "1_left",
    no_choice = 0) %>% 
  group_by(experiment, folder, choice_set) %>% 
  mutate(image = 1:n()) %>% 
  ungroup()


merged_right = merged_counts %>% 
  select(all_of(vars_to_select), ends_with("right")) %>% 
  set_names(c(vars_to_select, names_mapping$new_name)) %>% 
  mutate(
    alternative = "2_right",
    no_choice = 0) %>% 
  group_by(experiment, folder, choice_set) %>% 
  mutate(image = 1:n()) %>% 
  ungroup()


merged_no_choice =  merged_counts %>% 
  select(all_of(vars_to_select)) %>% 
  bind_cols(
    as.data.frame(matrix(0.0, nrow = nrow(merged_counts), ncol = length(names_mapping$new_name))) %>% 
      set_names(names_mapping$new_name)
  ) %>% 
  mutate(
    count = merged_counts$count_no_choice,
    alternative = "3_none",
    no_choice = 1) %>% 
  group_by(experiment, folder, choice_set) %>% 
  mutate(image = 1:n()) %>% 
  ungroup()


merged_all = merged_left %>% 
  bind_rows(merged_right) %>% 
  bind_rows(merged_no_choice) %>% 
  arrange(date_time, alternative) %>% 
  select(experiment, date_time, folder, choice_set, alternative, all_of(names_mapping$new_name), no_choice, image) %>% 
  mutate(alternative = as.factor(alternative)) %>% 
  filter(folder != "F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0")


merged_all %>% 
  filter(!complete.cases(.))


write_csv(merged_all, here("japanese_fly/out/counts_japanese_choice_format.csv"))


