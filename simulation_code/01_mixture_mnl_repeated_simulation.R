library(tidyverse)
library(mlogit)
library(opdesmixr)
library(here)

source(here("simulation_code/utils.R"))


q = 5
J = 2
S = 60


beta_vec_01 = rep(5, 14)


# Create an I-optimal design
design_array_i_opt_01_opdesmixr = opdesmixr::mnl_mixture_coord_exch(
  q = q, J = J, S = S, 
  beta = beta_vec_01, seed = 2021, order = 2, opt_crit = "I", transform_beta = F,
  n_random_starts = 1, n_cores = 1, max_it = 2)


# Repeat with 50 respondents
design_array_i_opt_01_df = lapply(1:50, function(r){
  print(r)
  
  S = dim(design_array_i_opt_01_opdesmixr$X)[3]
  
  out = design_array_i_opt_01_opdesmixr$X %>%
    mnl_design_array_to_dataframe() %>% 
    mutate(point = choice_set) %>% 
    mutate(
      choice_set = as.character(as.numeric(choice_set) + S*(r - 1)),
      resp = r
    )
  
  return(out)
  
}) %>% 
  bind_rows()


design_array_i_opt_01 = mnl_design_dataframe_to_array(design_array_i_opt_01_df %>% select(-resp, -point))



simulated_data_i_opt_01 = simulate_mixture_choice_data(design_array_i_opt_01, beta = beta_vec_01, order = 2, append_model_matrix = F)



# Transform data to be used in mlogit
sim_data_mlogit_i_opt_01 <- dfidx(
  simulated_data_i_opt_01 %>% 
    select(starts_with("comp"), choice_set, choice), 
  idx = c("choice_set"))

model_i_opt_01 = mlogit(choice ~ -1 + 
                          comp_1 + comp_2 + comp_3 + comp_4 + 
                          comp_1:comp_2 + comp_1:comp_3 + comp_1:comp_4 + comp_1:comp_5 +
                          comp_2:comp_3 + comp_2:comp_4 + comp_2:comp_5 +
                          comp_3:comp_4 + comp_3:comp_5 +
                          comp_4:comp_5,
                        sim_data_mlogit_i_opt_01)
model_i_opt_01


# My implementation of maximum likelihood estimation for MNL
my_ml_est_01 = maximum_likelihood_est(
  design_array = design_array_i_opt_01, 
  y = sim_data_mlogit_i_opt_01$choice, 
  order = 2, 
  starting_par = rep(0, length(beta_vec_01)), 
  maxit = 100)

my_ml_est_01$par

# Same results as mlogit
sum((my_ml_est_01$par - model_i_opt_01$coefficients)^2)




# Summarize data into counts per respondent
simulated_data_i_opt_01_counts = simulated_data_i_opt_01 %>% 
  mutate(resp = design_array_i_opt_01_df$resp,
         point = design_array_i_opt_01_df$point) %>% 
  group_by_at(vars(starts_with("comp"), point)) %>% 
  summarize(N_js = sum(choice)) %>% 
  mutate(point = as.numeric(point)) %>% 
  arrange(point) %>% 
  ungroup()


# Convert design into an array (the ML estimation function needs it)
simulated_data_i_opt_01_counts_array = simulated_data_i_opt_01_counts %>% 
  select(-N_js) %>% 
  rename(choice_set = point) %>% 
  mnl_design_dataframe_to_array()



# My implementation of maximum likelihood estimation for MNL with counts
my_ml_est_repeated_respondents_01 = maximum_likelihood_est_repeated_respondents(
  design_array = simulated_data_i_opt_01_counts_array, 
  Ny = simulated_data_i_opt_01_counts$N_js, 
  order = 2, 
  starting_par = rep(0, length(beta_vec_01)), 
  maxit = 1000)


my_ml_est_repeated_respondents_01$par


# Same results as mlogit
sum((my_ml_est_repeated_respondents_01$par - model_i_opt_01$coefficients)^2)

# And exactly the same result as with the "long" format (i.e., without the counts)
sum((my_ml_est_repeated_respondents_01$par - my_ml_est_01$coefficients)^2)



