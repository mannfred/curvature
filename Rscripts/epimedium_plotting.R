library(tidyverse)
library(here)

data_tbl <- 
    left_join(
        read.csv(
                here("data/epimedium_curv_curvature.csv"),
                header=TRUE) %>% 
        as.tibble(),
        read.csv(
                here("data/epimedium_curv_size_data.csv"),
                header=TRUE) %>%
        as.tibble() %>%
        slice(., 28:58), #isolate rows w E. koreanum data
        by="species_individual_panicle_flower" 
              ) %>%
    dplyr::select(4, 2, 3, 5, 6) #isolate columns 2:5 and rearrange


qplot(data=data_tbl, x=sepal_size_mm, y=adjusted_curvature) +
  ggtitle('curvature vs size')
   
        
      
