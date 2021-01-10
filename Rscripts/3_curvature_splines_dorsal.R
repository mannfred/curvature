library(curvr)
library(geomorph)
library(here)
library(Momit)
library(Momocs)
library(pracma)
library(tidyverse)

# -------------------------
# import complete landmark data


# automatically scaled by readmulti.tps()
epi_lmk_data <- 
  geomorph::readmulti.tps(
    c(here("data/epimedium_photos/koreanum/koreanum_appended_geomorph.TPS"),
      here("data/epimedium_photos/violaceum/violaceum_appended_geomorph.TPS")),
    specID="imageID")[,,-31] #remove K_1_2_1 (no size data)




# subset to dorsal landmarks (LMs 1-15, LM31 is the apex of the nectar spur to be omitted)
# Momocs::Ldk() converts to Coo object
dorsal <- epi_lmk_data[1:15,,]  %>%  Ldk()



# ------------
# alignment

# align the longest axis of each shape along the x-axis
dorsal_x <- coo_alignxax(dorsal)


# ---------------------
# calculate curvature


# extract the lower and upper x boundaries
baselines_list <- 
  dorsal_x$coo %>% 
  lapply(., function(b) c(unlist(b[,1])[1], unlist(b[,1])[15]))

# fit smoothing splines
spline_list <-
  dorsal_x$coo %>% 
  lapply(smooth.spline)

# compute first deriv of splines
deriv_list <-
  spline_list %>% 
  lapply(predict, deriv = 1)

# fit linear models to first derivs
model_list <- list()

for (i in 1:length(deriv_list)) {
  model_list[[i]] <- 
    lm(deriv_list[[i]][[2]] ~ deriv_list[[i]][[1]])
}


# compute curvature 
curvature_tbl <-
  mapply(spline_curvature, model_list, baselines_list) %>% 
  enframe() %>% 
  mutate(total_K = value*-180/pi) %>% 
  select(3)
  

# estimate arclength as perimeter
perim_list <- coo_perim(Ldk(dorsal_x))

# combine curv values with perimeter values
alltogether_tbl <- 
  perim_list %>% 
  unlist() %>% 
  enframe() %>% 
  mutate(total_K = curvature_tbl$total_K) %>% 
  mutate(name = dimnames(epi_lmk_data)[[3]]) %>% 
  rename(perimeter = value) 

# write rds
write_rds(alltogether_tbl, path = here("data/RDS_files/spline_curvature_tbl_dorsal.rds"))
