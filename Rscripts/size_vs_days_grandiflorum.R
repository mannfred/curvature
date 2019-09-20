library(tidyverse)
library(here)

#take "stages_days_sort" object from seq_alignment_*.R files and add size column

# # A tibble: 828 x 3
# ID                 elapsed_days stage
# <fct>                     <dbl> <chr>
#   1 grandiflorum_1_3_7            1 -    
#   2 grandiflorum_1_3_7            2 -    
#   3 grandiflorum_1_3_7            3 -    
#   4 grandiflorum_1_3_7            4 -    
#   5 grandiflorum_1_3_7            5 -    
#   6 grandiflorum_1_3_7            6 -    
#   7 grandiflorum_1_3_7            7 -    
#   8 grandiflorum_1_3_7            8 G    
#   9 grandiflorum_1_3_7            9 T    
#  10 grandiflorum_1_3_7           10 T    

#then plot elapsed_days vs size

stages_days_sort <- readRDS("stages_days_sort.rds")
