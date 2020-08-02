# measure curvature from epimedium photos

devtools::install("C:/Users/mannfred/Google Drive/UBC Botany/curvr", build=TRUE, quick=TRUE)

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

# subset to ventral landmarks (LMs 16-30, LM31 is the apex of the nectar spur to be omitted)
# Momocs::Ldk() converts to Coo object
ventral <- epi_lmk_data[16:30,,] %>%  Ldk()



# ------------
# alignment

# align the longest axis of each shape along the x-axis
ventral_x <- coo_alignxax(ventral)



# --------------------------

# compute circularity for each sample
# ~ 1 = circle, higher numbers = non-circles
circularity <-  1/as.numeric(coo_circularitynorm(ventral_x))

# saveRDS(circularity, file=here('data/RDS_files/circularity.rds'))



# -------------------------
# investigate polynomial model fits

# examine R2 for different polynomial models
# R2 for degree 1 to 6
r <- numeric()
for (i in 1:10) { r[i] <- npoly(ventral_p$coords[,,7], degree=i)$r2 }
plot(2:10, r[2:10], type='b', pch=20, col='red', main='R2 / degree')


# try out some polynomials
poly_list <- 
  map(ventral_x$coo, Momocs::npoly, degree = 4)



# plot landmarks
coo_plot(ventral_x$coo[32])

# inspect polynomial fit
# inspect 4, 13, 14, 21, 22, 23..
poly_list[[32]] %>% 
  npoly_i() %>% #calculates shape from polynomial model
  coo_draw(border='red')


# extract all R2 fits
r2 <- list()

for (i in 1:length(poly_list)) {
  r2[[i]] <- poly_list[[i]]$r2 
}

# visualize R2 distribution
# lowest R2 is 0.97 for degrees = 4
plot(1:57, as.numeric(r2), title(main='degrees=4'))


# ---------------------
# calculate curvature


# extract the lower and upper bounds from b[4] and b[5], respectively
baselines_list <- 
  poly_list %>% 
  lapply(., function(b) c(unlist(b[4])[1], unlist(b[5])[1]))


# calculate total curvature 
# units are in degrees per unit length
curvature_tbl <- 
  curvr::total_curvature %>% 
  mapply(., baselines_list, poly_list, 500) %>% 
  enframe()

saveRDS(curvature_tbl, file = here("data/RDS_files/curvature_tbl_ventral.rds"))



# ---------------------------
# compute arc length
# either from polynomials or from shapes' perimeter

# a list of min bounds, '[[' is a subsetting function
min_baselines <- 
  baselines_list %>% 
  sapply(., "[[", 1)

# a list of max bounds
max_baselines <- 
  baselines_list %>% 
  sapply(., "[[", 2)

# calculate arc length by first parameterizing polynomials by t
param_list <- 
  lapply(poly_list, curvr::parameterize) 

# calculate arclength for every bounded polynomial
lengths_list <- 
  mapply(pracma::arclength, param_list, min_baselines, max_baselines) %>% 
  as_tibble() %>% 
  slice(., 1) %>% # keep only row 1
  as.list()


# or estimate arclength as perimeter
perim_list <- coo_perim(Ldk(ventral))


# divide curvature by arclength 
alltogether_tbl <- 
  perim_list %>% 
  unlist() %>% 
  enframe() %>% 
  left_join(., curvature_tbl, by = 'name') %>%  # join arclength columns and curvature columns 
  rename(arclength = value.x, total_curvature = value.y) %>% 
  mutate(adjusted_curvature = total_curvature/arclength) %>%  # new column is adjusted curvature
  mutate(name = dimnames(epi_lmk_data)[[3]])

write.csv(alltogether_tbl, here('data/epimedium_adj_curvature_ventral.csv')) 

