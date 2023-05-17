library(tidyverse)
library(here)

# With the shitty One Drive app I can have this locally
# experiment_folder = "/Volumes/hdd_1tb/OneDrive - KU Leuven/THESIS STATISTICS/EXPERIMENTS/"
experiment_folder = "C:/Users/SET-L-DI-02608/OneDrive - KU Leuven/THESIS STATISTICS/EXPERIMENTS/"

files_folders_all = grep("design", list.files(experiment_folder, full.names = T), value = T, ignore.case = T)

files_folders = grep("japanese_[0-9][0-9]", files_folders_all, value = T, perl = T)


all_counts = lapply(1:length(files_folders), function(i){
  subfolders_i = grep("202", list.files(files_folders[i], full.names = F), value = T)
  
  print(files_folders[i])
  
  design_filename = grep("design", list.files(files_folders[i]), value = T)
  
  design_exp_i = read.table(paste0(files_folders[i], "/", design_filename), header = T)
  
  key_df_i = design_exp_i %>% 
    select(5:ncol(.)) %>% 
    mutate_at(vars(starts_with("Var")), ~as.character(round(.))) %>% 
    mutate_at(vars(!starts_with("Var")), ~substr(as.character(.), 1, 1))
  
  keys = apply(key_df_i, 1, function(x) paste0(x, collapse = ""))
  
  names_colors = design_exp_i %>% 
    select(5:ncol(.)) %>% 
    select(!starts_with("Var")) %>% 
    names()
  
  design_exp_i_2 = design_exp_i %>% 
    select(starts_with("Var")) %>% 
    set_names(c(paste0(names_colors[1:(length(names_colors)/2)], "_left"), 
                paste0(names_colors[1:(length(names_colors)/2)], "_right"))) %>% 
    mutate(key = keys) %>% 
    distinct()
  
  out_exp_i = lapply(1:length(subfolders_i), function(j){
    subfolder_ij = list.files(paste0(files_folders, "/", subfolders_i[j]), full.names = T)[1]
    date_ij = subfolders_i[j]
    subfolders_data_ij = list.files(subfolder_ij, full.names = F)
    
    out_ij = lapply(1:length(subfolders_data_ij), function(k){
      counts_file_ijk = paste0(subfolder_ij, "/", subfolders_data_ij[k], "/Exp_140000/resulst.txt")
      
      if(file.exists(counts_file_ijk)){
        data = read.table(counts_file_ijk, header = F) %>% 
          mutate(folder = subfolders_data_ij[k]) %>% 
          as_tibble()
      } else {
        data = NULL
      }
      
      return(data)
    }) %>% 
      bind_rows() %>% 
      mutate(date = date_ij)
  }) %>% 
    bind_rows() %>% 
    mutate(experiment = i) %>% 
    left_join(design_exp_i_2, by = c("folder" = "key"))
  
}) %>% 
  bind_rows() %>% 
  rename(time_of_day = V1, count_left = V2, count_right = V3) %>% 
  arrange(date, time_of_day) %>% 
  mutate(ix = 1:nrow(.)) %>% 
  mutate(experiment = as.integer(fct_reorder(as.factor(experiment), ix))) %>%  # reorder experiment number by date
  select(-ix)





# Check if there are missing values
all_counts %>% filter(!complete.cases(.))


# Create the choice set variable
all_counts$choice_set = NA_integer_
all_counts$choice_set[1] = 1
counter = 1
for(i in 2:nrow(all_counts)){
  if(all_counts$folder[i] != all_counts$folder[i-1]){
    counter = counter + 1
  }
  all_counts$choice_set[i] = counter
}



date_char = as.character(Sys.Date())
write_csv(all_counts, here(paste0("data/counts_japanese_", date_char, ".csv")))

