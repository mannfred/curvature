# measure curvature from epimedium photos

install("C:/Users/mannfred/Google Drive/UBC Botany/curvr", build=TRUE, quick=TRUE)
library(curvr)
library(here)
library(Momit)
library(Momocs)
library(pracma)
library(tidyverse)

# why not import both kore and viol here?
# import data 
 dorsal_data <- 
   here('data/epimedium_photos/koreanum/koreanum_dorsal_appended.TPS') %>%
   import_tps()


# format data and apply scaling 
dorsal_lst <- 
  map2(dorsal_data$coo, dorsal_data$scale, function(x, y) x*y) %>% 
  Ldk()


# fit polynomials to landmarks
dorsalcurv_lst <- 
  map(dorsal_lst$coo, Momocs::npoly, 3)


# calculate arc length
paramfun_lst <- 
  lapply(dorsalcurv_lst, curvr::as_param)  # a list of 31 parameterized polynomial functions


# extract the lower and upper bounds from b[5] and b[4], respectively
baselines_lst <- 
  dorsalcurv_lst %>% 
  lapply(., function(b) c(unlist(b[5])[1], unlist(b[4])[1]))

# a list of min bounds, '[[' is a subsetting function
min_baselines <- 
  baselines_lst %>% 
  sapply(., "[[", 1)

# a list of max bounds
max_baselines <- 
  baselines_lst %>% 
  sapply(., "[[", 2)

# calculates arclength for every bounded polynomial
lengths_lst <- 
  mapply(arclength, paramfun_lst, min_baselines, max_baselines) %>% 
  as_tibble() %>% 
  slice(., 1) %>% # keep only row 1
  as.list()


# calculate total curvature 
# units are in degrees per unit length
curvature_tbl <- 
  curvr::totalK %>% 
  mapply(., baselines_lst, dorsalcurv_lst, 500) %>% 
  as_tibble() %>% 
  gather()

saveRDS(curvature_tbl, file = here("data/RDS_files/curvature_tbl.rds"))



# merge curvature and arclength information
alltogether_tbl <- # list of arclengths
  lengths_lst %>% 
  unlist() %>% 
  as.tibble() %>% # unpivot column names to row names
  gather() %>% 
  # bind arclength column and curvature column 
  # mask Momocs::select 
  # fetch ID tags from data.csv isolate ID column 20:50 for
  # koreanum, 51:77 for violaceum, 1:19 for grandiflorum
  {bind_cols(dplyr::select(., value), 
             dplyr::select(curvature_tbl, value), 
             dplyr::select(read.csv(here("data/epimedium_curv_size_data.csv"), header = TRUE) %>%
             dplyr::select(species_individual_panicle_flower) %>% 
             slice(., 20:50) %>% 
             as.tibble(), species_individual_panicle_flower))
  } %>% 
  rename(arclength = value, total_curvature = value1) %>% #old colnames were 'value' and 'value1'
  mutate(adjusted_curvature = total_curvature/arclength) #new column is adjusted curvature



 write.csv(., here('data/epimedium_curvature_koreanum.csv')) 







