library(curvr)
library(geomorph)
library(here)
library(Momit)
library(Momocs)
library(pracma)
library(splines2)
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
# align w x-axis

# align the longest axis of each shape along the x-axis
dorsal_x <- coo_alignxax(dorsal)



# ----------------------------
# plot fitted splines onto landmarks

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


# -----------------------------
# compute curvature of fitted splines
x1 <- seq(1, 10, 0.2)
y1 <- x1^2 
mat <- data.frame(x = x1, y =y1)

s0 <- smooth.spline(mat)

# first deriv values
s1 <- predict(s0, deriv=1)

# fit spline func to first deriv values
s1func <- splinefun(x=x1, y=s1$y)

# second deriv func
s2func <- splinefun(x=x1,  y = s1func(x1, deriv = 1)) 


k_fun <- function(x) {
  f1 <- s1func
  f2 <- s2func
  ((f2(x))/((1 + (f1(x)^2))^1.5)) * (sqrt(1 + (f1(x))^2))
}

integrate(k_fun, 1, 10)
