library(rstan)
library(tidyverse)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)

design_mapping = read_csv("data/design_mapping.csv")

# There's a mistake here, should have filtered out experiment number 1.
counts_european = read_csv(here("out/counts_european_choice_format.csv"))

max_counts = counts_european %>% 
  group_by(experiment, choice_set, folder, alternative) %>% 
  summarize(count = max(count)) %>% 
  ungroup %>% 
  filter(alternative != "3_none")




count_data_flies = max_counts %>% 
  full_join(
    counts_european %>% 
      select(experiment, alternative, choice_set, folder, R, O, Y, G, B, P, UV, intensity, no_choice) %>% 
      distinct()
  ) %>% 
  arrange(experiment, choice_set, folder, alternative) %>% 
  group_by(experiment, choice_set, folder) %>% 
  mutate(no_choice_count = 80 - sum(count, na.rm = T)) %>% 
  mutate(count = ifelse(is.na(count), no_choice_count, count)) %>% 
  select(-no_choice_count) %>% 
  ungroup()






unique_experiment_choice_set_01 = count_data_flies %>% 
  select(experiment, choice_set) %>% 
  distinct()

index_dataframe_01 = lapply(1:nrow(unique_experiment_choice_set_01), function(i){
  exp_i = unique_experiment_choice_set_01$experiment[i]
  cs_i =  unique_experiment_choice_set_01$choice_set[i]
  
  indices_i = which(count_data_flies$choice_set == cs_i & count_data_flies$experiment == exp_i)
  out = tibble(
    choice_set = cs_i, 
    start = indices_i[1], 
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>% 
  bind_rows()






X_main <- count_data_flies %>%
  select(c('R','O','Y','G','B','P','UV','intensity','no_choice')) %>%
  as.data.frame()

X_color_pairwise <- combn(c('R','O','Y','G','B','P','UV'), 2, function(x) X_main[x[1]]*X_main[x[2]], simplify = FALSE) %>%
  bind_cols() %>%
  set_names(c('R*O','R*Y','R*G','R*B','R*P','R*UV',
              'O*Y','O*G','O*B','O*P','O*UV',
              'Y*G','Y*B','Y*P','Y*UV',
              'G*B','G*P','G*UV',
              'B*P','B*UV',
              'P*UV'))

X_color_intensity <- X_main %>%
  select(c('R','O','Y','G','B','P','UV','intensity')) %>%
  mutate(across(c('R','O','Y','G','B','P','UV','intensity'), ~.*intensity)) %>%
  set_names(c('R*intensity','O*intensity','Y*intensity','G*intensity',
              'B*intensity','P*intensity','UV*intensity','intensity^2'))






X_stan_01 <- X_main %>% 
  select(-UV, -intensity) %>% 
  bind_cols(X_color_pairwise) %>% 
  bind_cols(X_color_intensity) %>% 
  as.matrix()





stan_data_01 <- list(
  J = nrow(count_data_flies)/max(count_data_flies$choice_set), 
  S = max(count_data_flies$choice_set), 
  n_var = ncol(X_stan_01), 
  n_obs = nrow(X_stan_01),
  Ny = count_data_flies$count,
  X = X_stan_01,
  start = index_dataframe_01$start, # the starting observation for each decision/choice set
  end = index_dataframe_01$end  # the ending observation for each decision/choice set
)


# index_dataframe = lapply(unique(simulated_data_counts$choice_set), function(i){
#   indices_i = which(simulated_data_counts$choice_set == i)
#   out = tibble(
#     choice_set = i, start = indices_i[1], end = indices_i[length(indices_i)]
#   )
# }) %>% 
#   bind_rows()
# 
# 
# stan_data_02 <- list(
#   J = J+1, 
#   S = max(index_dataframe$choice_set), 
#   n_var = ncol(X_stan_02), 
#   n_obs = nrow(X_stan_02),
#   Ny = simulated_data_counts$Ny,
#   X = X_stan_02,
#   start = index_dataframe$start, # the starting observation for each decision/choice set
#   end = index_dataframe$end  # the ending observation for each decision/choice set
# )


m_01 = "
data {
    int<lower=2> J; // number of alternatives/outcomes
    int<lower=2> S; // number of choice sets
    int<lower=1> n_var; // of covariates
    int<lower=1> n_obs; // of observations
    vector[n_obs] Ny; // counts of the response variable
    matrix[n_obs, n_var] X; // attribute matrix
    int start[S]; // the starting observation for each decision/choice set
    int end[S]; // the ending observation for each decision/choice set
}

parameters {
    vector[n_var] beta;  // attribute effects 
}

model {
    
    for(s in 1:S){
    
      vector[J] utilities;
      
      utilities = X[start[s]:end[s]]*beta;
      
      target += sum(log_softmax(utilities) .* Ny[start[s]:end[s]]);
    }
  
}
"

# 1 minute for 2000 iterations and 4 chains in 4 cores
model_stan_01 <- stan(
  model_code = m_01, 
  data = stan_data_01, 
  warmup = 2000,
  iter = 5000, 
  chains = 4, 
  seed = 2023,
  cores = 4)

model_stan_01




summary(model_stan_01, pars = c("beta"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_01)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



summary(model_stan_01, pars = c("beta"), probs = c(0.5))$summary %>% 
  as.data.frame() %>% 
  mutate(variable = colnames(X_stan_01)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), color = "dark grey") +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



# To use as prior in design
summary(model_stan_01, pars = c("beta"), probs = c(0.5))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "50%") %>% 
  mutate(variable = colnames(X_stan_01))















count_data_flies_wo_nochoice = count_data_flies %>% 
  filter(no_choice != 1)

unique_experiment_choice_set_02 = count_data_flies_wo_nochoice %>% 
  select(experiment, choice_set) %>% 
  distinct()

index_dataframe_02 = lapply(1:nrow(unique_experiment_choice_set_02), function(i){
  exp_i = unique_experiment_choice_set_02$experiment[i]
  cs_i =  unique_experiment_choice_set_02$choice_set[i]
  
  indices_i = which(count_data_flies_wo_nochoice$choice_set == cs_i & count_data_flies_wo_nochoice$experiment == exp_i)
  out = tibble(
    choice_set = cs_i, 
    start = indices_i[1], 
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>% 
  bind_rows()






X_main <- count_data_flies_wo_nochoice %>%
  select(c('R','O','Y','G','B','P','UV','intensity')) %>%
  as.data.frame()

X_color_pairwise <- combn(c('R','O','Y','G','B','P','UV'), 2, function(x) X_main[x[1]]*X_main[x[2]], simplify = FALSE) %>%
  bind_cols() %>%
  set_names(c('R*O','R*Y','R*G','R*B','R*P','R*UV',
              'O*Y','O*G','O*B','O*P','O*UV',
              'Y*G','Y*B','Y*P','Y*UV',
              'G*B','G*P','G*UV',
              'B*P','B*UV',
              'P*UV'))

X_color_intensity <- X_main %>%
  select(c('R','O','Y','G','B','P','UV','intensity')) %>%
  mutate(across(c('R','O','Y','G','B','P','UV','intensity'), ~.*intensity)) %>%
  set_names(c('R*intensity','O*intensity','Y*intensity','G*intensity',
              'B*intensity','P*intensity','UV*intensity','intensity^2'))






X_stan_02 <- X_main %>% 
  select(-UV, -intensity) %>% 
  bind_cols(X_color_pairwise) %>% 
  bind_cols(X_color_intensity) %>% 
  as.matrix()





stan_data_02 <- list(
  J = nrow(count_data_flies_wo_nochoice)/max(count_data_flies_wo_nochoice$choice_set), 
  S = max(count_data_flies_wo_nochoice$choice_set), 
  n_var = ncol(X_stan_02), 
  n_obs = nrow(X_stan_02),
  Ny = count_data_flies_wo_nochoice$count,
  X = X_stan_02,
  start = index_dataframe_02$start, # the starting observation for each decision/choice set
  end = index_dataframe_02$end  # the ending observation for each decision/choice set
)


# 1 minute for 2000 iterations and 4 chains in 4 cores
model_stan_02 <- stan(
  model_code = m_01, 
  data = stan_data_02, 
  warmup = 2000,
  iter = 5000, 
  chains = 4, 
  seed = 2023,
  cores = 4)

model_stan_02




summary(model_stan_02, pars = c("beta"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_02)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")


summary(model_stan_02, pars = c("beta"), probs = c(0.5))$summary %>% 
  as.data.frame() %>% 
  mutate(variable = colnames(X_stan_02)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), color = "dark grey") +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



# To use as prior in design
summary(model_stan_02, pars = c("beta"), probs = c(0.5))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "50%") %>% 
  mutate(variable = colnames(X_stan_02))

