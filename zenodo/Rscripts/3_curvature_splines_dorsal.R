# Fit splines and calculate total curvature from landmarked specimens


# devtools::install_github("mannfred/curvr")
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
    c(here("data/raw_data/epimedium_photos/koreanum/koreanum_appended_geomorph.TPS"),
      here("data/raw_data/epimedium_photos/violaceum/violaceum_appended_geomorph.TPS")),
    specID="imageID")[,,-31] #remove K_1_2_1 (no size data)



# subset to dorsal landmarks (LMs 1-15, LM31 is the apex of the nectar spur)
# Momocs::Ldk() converts to Coo object
dorsal <- epi_lmk_data[c(31,1:15),,]  %>%  Ldk()

# align the longest axis of each shape along the x-axis
dorsal_x <- coo_alignxax(dorsal)

# simplify nested coo objects
dorsal_xsimple <- unlist(dorsal_x, recursive=F) 



# ----------------------------------
# calculate total curvature by curvr

# split the curves in two
one <- lapply(dorsal_xsimple, function(x) x[][1:7,])
two <- lapply(dorsal_xsimple, function(x) x[][8:16,])

# realign with x-axis
onex <- coo_alignxax(Ldk(one)) %>% unlist(., recursive = F)
twox <- coo_alignxax(Ldk(two)) %>% unlist(., recursive = F)


# extract the lower and upper x boundaries
baselines_onex <- 
  onex %>% 
  lapply(., function(b) c(unlist(b[,1])[1], unlist(b[,1])[7]))

baselines_twox <- 
  twox %>% 
  lapply(., function(b) c(unlist(b[,1])[1], unlist(b[,1])[9]))


# fit interpolating splines and compute total curvature using package curvr
curvature_onex <-
 mapply(curvature_spline, onex, baselines_onex, 'ip') %>% 
  enframe() %>% 
  mutate(total_K = abs(value)*(180/pi)) 

curvature_twox <-
  mapply(curvature_spline, twox, baselines_twox, 'ip') %>% 
  enframe() %>% 
  mutate(total_K = abs(value)*(180/pi)) 


# compile data into one table
curv_data <- 
  tibble(
    total_K = curvature_onex$total_K + curvature_twox$total_K,
    name = dimnames(epi_lmk_data)[[3]]) 

# write rds
# write_rds(curv_data, path = here("data/derived_data/RDS_files/spline_curvature_tbl_dorsal.rds"))
# to compare to other metrics in "7_compare_altmetrics.R"

# ----------------------------------
# plotting splines
coords <- matrix()
s0 <- vector()

for (i in 1:length(onex)){
  # test coords
  coords <- onex[[i]]
  
  s0 <- spline(coords)
  
  # plot
  plot(coords)
  lines(s0, col='red', lw=2)
  title(main = paste('specimen #', i, sep=""))
}

