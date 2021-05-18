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
dorsal <- epi_lmk_data[1:15,,]  %>%  Ldk()

# align the longest axis of each shape along the x-axis
dorsal_x <- coo_alignxax(dorsal)


# -------------------------
# fit circles to landmarks

# estimate radius of circle
# extract column 1 from a list of matrices
xcoords <- lapply(dorsal_x$coo, function(i) i[,1])

# extract column 2 from a list of matrices
ycoords <- lapply(dorsal_x$coo, function(i) i[,2])

# fit circles 
# RMS error is printed to console, 
# so it has been saved to circle_rms_error.csv
circles <- mapply(pracma::circlefit, xcoords, ycoords)

# extract radii
radius <- circles[3,]

# save to .csv 
# write.csv(radius, file=here('data/derived_data/circles_fitby_pracma.csv'))
# for use in "7_compare_altmetrics.R"
# NOTE: rms error is directly printed onto console by circlefit()
# but not stored in the object "circles". It's been pasted into "circle_rms_error_pracma.csv".

# visualize circle fits
coo_plot(dorsal_x$coo[4])

x0 <- circles[1,4]
y0 <- circles[2,4]
r0 <- circles[3,4]

w  <- seq(0, 2*pi, len=100)
xx <- r0 * cos(w) + x0
yy <- r0 * sin(w) + y0
lines(xx, yy, col="red")


# --------------------------------
# calculate curvature


# extract the lower and upper x boundaries
baselines_list <- 
  dorsal_x$coo %>% 
  lapply(., function(b) c(unlist(b[,1])[1], unlist(b[,1])[15]))

# for nested coo objects
dorsal_xsimple <- unlist(dorsal_x, recursive=F) 


# fit splines and compute curvature
curvature_tbl <-
 mapply(curvature_spline, dorsal_xsimple, baselines_list) %>% 
  enframe() %>% 
  mutate(total_K = abs(value)*(180/pi)) 


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
# write_rds(alltogether_tbl, path = here("data/derived_data/RDS_files/spline_curvature_tbl_dorsal.rds"))



# plotting splines
coords <- matrix()
s1 <- vector()

for (i in 1:length(dorsal_x)){
  # test coords
  coords <- dorsal_x$coo[[i]]
  
  s1 <- smooth.spline(coords)
  
  # plot
  plot(coords)
  lines(s1)
}

