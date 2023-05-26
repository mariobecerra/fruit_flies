library(tidyverse)
library(here)

counts_european = read_csv(here("out/counts_european_choice_format.csv"))


## Experiment 5

experiment_5_left_vs_right = counts_european %>% 
  filter(experiment == 5) %>% 
  filter(alternative != "3_none") %>% 
  mutate(UV_on = ifelse(UV > 0, "UV on", "UV off")) %>% 
  select(date_time, choice_set, alternative, count, UV_on) %>% 
  pivot_wider(names_from = alternative, values_from = count) %>% 
  rename(
    left = `1_left`,
    right = `2_right`) %>% 
  mutate(right_minus_left = right - left)

mean(experiment_5_left_vs_right$right_minus_left)
median(experiment_5_left_vs_right$right_minus_left)

experiment_5_left_vs_right %>% 
  ggplot() +
  geom_histogram(aes(right_minus_left)) +
  geom_vline(xintercept = 0) +
  theme_bw()



experiment_5_left_vs_right %>% 
  ggplot() +
  geom_histogram(aes(right_minus_left)) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  facet_wrap(~UV_on, scales = "free")


experiment_5_left_vs_right %>% 
  group_by(UV_on) %>% 
  summarize(
    mean_right_minus_left = mean(right_minus_left),
    median_right_minus_left = median(right_minus_left)
  )


experiment_5_left_vs_right %>% 
  mutate(ix = 1:nrow(.)) %>% 
  ggplot() +
  geom_point(aes(x = right_minus_left, y = ix, color = right_minus_left), size = 0.2) +
  geom_segment(aes(x = 0, xend = right_minus_left, y = ix, yend = ix, color = right_minus_left)) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, nrow(experiment_5_left_vs_right), by = 100)) +
  geom_hline(yintercept = 716)


# # Not sure if this plot says a lot
# experiment_5_left_vs_right %>%
#   # head(500) %>%
#   mutate(ix = 1:nrow(.)) %>%
#   ggplot() +
#   geom_point(aes(x = left, y = ix, color = right_minus_left), size = 0.2) +
#   geom_point(aes(x = right, y = ix, color = right_minus_left), size = 0.2) +
#   geom_segment(aes(x = left, xend = right, y = ix, yend = ix, color = right_minus_left)) +
#   theme_bw()




colors = c("R", "O", "Y", "G", "B", "P")


aaa = counts_european %>% 
  filter(experiment == 5) %>% 
  filter(alternative != "3_none") %>% 
  select(date_time, all_of(colors)) %>% 
  distinct() %>% 
  select(all_of(colors))

color_vector = rep(NA_character_, nrow(aaa))
for(i in seq_along(color_vector)){
  color_char = names(aaa)[which(as.numeric(aaa[i,]) > 0)]
  if(length(color_char) == 1) color_vector[i] = color_char
}





experiment_5_left_vs_right_2 = experiment_5_left_vs_right %>% 
  mutate(ix = 1:nrow(.)) %>% 
  mutate(color = ifelse(is.na(color_vector), "UV", color_vector)) %>% 
  group_by(choice_set) %>% 
  mutate(init_ix_choice_set = min(ix),
         end_ix_choice_set = max(ix),
         center_ix_choice_set = median(ix)) %>% 
  ungroup()


experiment_5_left_vs_right_text =  experiment_5_left_vs_right %>% 
  mutate(ix = 1:nrow(.)) %>% 
  mutate(color = ifelse(is.na(color_vector), "UV", color_vector)) %>% 
  group_by(choice_set, color) %>% 
  summarize(init_ix_choice_set = min(ix),
         end_ix_choice_set = max(ix),
         center_ix_choice_set = median(ix)) %>% 
  ungroup()



experiment_5_left_vs_right_2 %>% 
  mutate(ix = 1:nrow(.)) %>% 
  ggplot() +
  geom_point(aes(x = right_minus_left, y = ix, color = right_minus_left), size = 0.2) +
  geom_segment(aes(x = 0, xend = right_minus_left, y = ix, yend = ix, color = right_minus_left)) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, nrow(experiment_5_left_vs_right), by = 100)) +
  geom_hline(aes(yintercept = init_ix_choice_set),
             data = experiment_5_left_vs_right_text) +
  geom_text(x = 20, 
            aes(
              y = center_ix_choice_set, 
              label = color
            ),
            data = experiment_5_left_vs_right_text,
            inherit.aes = F)









#### Experiments 1 to 4

experiments_1_to_4_left_vs_right = counts_european %>% 
  filter(experiment != 5) %>% 
  filter(alternative != "3_none") %>% 
  select(experiment, date_time, choice_set, alternative, count) %>% 
  pivot_wider(names_from = alternative, values_from = count) %>% 
  rename(
    left = `1_left`,
    right = `2_right`) %>% 
  mutate(right_minus_left = right - left)



mean(experiments_1_to_4_left_vs_right$right_minus_left)
median(experiments_1_to_4_left_vs_right$right_minus_left)

experiments_1_to_4_left_vs_right %>% 
  ggplot() +
  geom_histogram(aes(right_minus_left)) +
  geom_vline(xintercept = 0) +
  theme_bw()



experiments_1_to_4_left_vs_right %>% 
  ggplot() +
  geom_histogram(aes(right_minus_left)) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  facet_wrap(~experiment, scales = "free")



experiments_1_to_4_left_vs_right %>% 
  group_by(experiment) %>% 
  summarize(
    mean_right_minus_left = mean(right_minus_left),
    median_right_minus_left = median(right_minus_left)
  )





experiments_1_to_4_left_vs_right %>% 
  mutate(ix = 1:nrow(.)) %>% 
  ggplot() +
  geom_point(aes(x = right_minus_left, y = ix, color = right_minus_left), size = 0.2) +
  geom_segment(aes(x = 0, xend = right_minus_left, y = ix, yend = ix, color = right_minus_left)) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, nrow(experiments_1_to_4_left_vs_right), by = 500)) 




experiments_1_to_4_left_vs_right_text =  experiments_1_to_4_left_vs_right %>% 
  mutate(ix = 1:nrow(.)) %>% 
  group_by(experiment) %>% 
  summarize(init_ix_choice_set = min(ix),
            end_ix_choice_set = max(ix),
            center_ix_choice_set = median(ix)) %>% 
  ungroup()



experiments_1_to_4_left_vs_right %>% 
  mutate(ix = 1:nrow(.)) %>% 
  ggplot() +
  geom_point(aes(x = right_minus_left, y = ix, color = right_minus_left), size = 0.2) +
  geom_segment(aes(x = 0, xend = right_minus_left, y = ix, yend = ix, color = right_minus_left)) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, nrow(experiments_1_to_4_left_vs_right), by = 500)) +
  geom_hline(aes(yintercept = init_ix_choice_set),
             data = experiments_1_to_4_left_vs_right_text) +
  geom_text(x = 15, 
            aes(
              y = center_ix_choice_set, 
              label = paste0("Experiment ", experiment)
            ),
            data = experiments_1_to_4_left_vs_right_text,
            inherit.aes = F)










