# measure curvature from epimedium photos

install("C:/Users/mannfred/Google Drive/UBC Botany/curvr", build=TRUE, quick=TRUE)

library(curvr)
library(here)
library(Momit)
library(Momocs)
library(pracma)
library(tidyverse)


# import data 
 kore_data <- 
   here('data/epimedium_photos/koreanum/koreanum_dorsal_appended.TPS') %>%
   import_tps()
 
 viol_data <-
   here('data/epimedium_photos/violaceum/violaceum_dorsal_appended.TPS') %>% 
   import_tps()

# combine
dorsal_data <-
  map2(kore_data, viol_data, c)


# scale
scaled <-
   map2(dorsal_data$coo, dorsal_data$scale, function(x, y) x*y) %>% 
   Ldk()
 
 
# fit polynomials to landmarks
poly_list <- 
 map(scaled$coo, Momocs::npoly, 3)


# calculate total curvature 
# units are in degrees per unit length
curvature_tbl <- 
  curvr::totalK %>% 
  mapply(., baselines_list, poly_list, 500) %>% 
  enframe()

saveRDS(curvature_tbl, file = here("data/RDS_files/curvature_tbl.rds"))


# calculate arc length by first parameterizing polynomials by t
param_list <- 
  lapply(poly_list, curvr::as_param)  # a list of 31 parameterized polynomial functions


# extract the lower and upper bounds from b[5] and b[4], respectively
baselines_list <- 
  poly_list %>% 
  lapply(., function(b) c(unlist(b[5])[1], unlist(b[4])[1]))

# a list of min bounds, '[[' is a subsetting function
min_baselines <- 
  baselines_list %>% 
  sapply(., "[[", 1)

# a list of max bounds
max_baselines <- 
  baselines_list %>% 
  sapply(., "[[", 2)

# calculate arclength for every bounded polynomial
lengths_list <- 
  mapply(arclength, param_list, min_baselines, max_baselines) %>% 
  as_tibble() %>% 
  slice(., 1) %>% # keep only row 1
  as.list()


# divide curvature by arclength 
alltogether_tbl <- 
  lengths_list %>% 
  unlist() %>% 
  enframe() %>% 
  left_join(., curvature_tbl, by = 'name') %>%  # join arclength columns and curvature columns 
  rename(arclength = value.x, total_curvature = value.y) %>% 
  mutate(adjusted_curvature = total_curvature/arclength) # new column is adjusted curvature


write.csv(alltogether_tbl, here('data/epimedium_adj_curvature.csv')) 

