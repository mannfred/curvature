library(tidyverse)
library(here)

data_tbl <- 
    left_join(
        read.csv(
                here("data/epimedium_curvature_violaceum.csv"),
                header=TRUE) %>% 
            as.tibble(),
        read.csv(
                here("data/epimedium_curv_size_data.csv"),
                header=TRUE) %>%
            as.tibble() %>%
        slice(., 51:77), #20:50 for koreanum, 51:77 for violaceum, 1:21 for grandiflorum
        by="species_individual_panicle_flower" 
              ) %>%
    dplyr::select(4, 2, 3, 5, 6) #isolate columns 2:5 and rearrange


qplot(data=data_tbl, x=sepal_size_mm, y=adjusted_curvature) +
  ggtitle('curvature vs size')

   
        
      
