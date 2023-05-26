# Creates a design with all six colors to check whether they attract roughly the same number of flies.
# With and without UV light.
# Randomized order of colors.

library(dplyr)
library(here)

set.seed(2023)

colors = c("R", "O", "Y", "G", "B", "P")
n_colors = length(colors)

random_order_1 = c(sample(colors, n_colors))
random_order_2 = c(sample(colors, n_colors))

df_aux = as.data.frame(matrix(0.0, ncol = n_colors + 1, nrow = 2*n_colors + 1))

names(df_aux) = c(colors, "UV")

df_aux$UV[(n_colors+1):nrow(df_aux)] = 1

for(i in 1:n_colors){
  df_aux[i, random_order_1[i]] = 1
  df_aux[i + n_colors, random_order_2[i]] = 1
}


design_df = df_aux %>% 
  mutate(choice_set = 1:nrow(df_aux)) %>% 
  bind_rows(
    df_aux %>% 
      mutate(choice_set = 1:nrow(df_aux))) %>% 
  arrange(choice_set)


write.csv(design_df, here("out/design_same_colors_UV.csv"), col.names = F)

