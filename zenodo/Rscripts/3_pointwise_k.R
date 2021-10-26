# Estimate point-curvature at each of the 16 landmarks
# and compare to the constant K value estimated from the 1/R metric

library(geomorph)
library(here)
library(Momit)
library(Momocs)
library(pracma)
library(tidyverse)


# ------------------------------------------
# import landmark data


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



# ----------------------------------------------------
# import circledata from "3_circle_fitting.R"

circledata <- read.csv(file=here('data/derived_data/circles_fitby_pracma.csv'))
arcs <- circledata$arcs
radius <- circledata$radius
csize <- circledata$csize



# ---------------------------------------------------
# calculate and plot point-wise K for each specimen

# split the curves in two
one <- lapply(dorsal_xsimple, function(x) x[][1:7,])
two <- lapply(dorsal_xsimple, function(x) x[][8:16,])

# realign with x-axis
onex <- coo_alignxax(Ldk(one)) %>% unlist(., recursive = F)
twox <- coo_alignxax(Ldk(two)) %>% unlist(., recursive = F)



ydis <- numeric(); 

for (i in 1:length(dorsal_x)){

# first curve segment
s0 <- spline(onex[[i]])

# create a spline function from coordinates
s0func <- splinefun(s0$x, s0$y)

# estimate y coords of first derivative
s1 <- s0func(s0$x, deriv = 1)

# create a function for first derivative
s1func <- splinefun(s0$x, s1)

# create a function for second derivative
s2func <- splinefun(x = s0$x, y = s1func(s0$x, deriv = 1))

# estimate points from second deriv function
s2coords <- s2func(onex[[i]][,1])

# second curve segment
s3 <- spline(twox[[i]])
s3func <- splinefun(s3$x, s3$y)
s4 <- s3func(s3$x, deriv = 1)
s4func <- splinefun(s3$x, s4)
s5func <- splinefun(x = s3$x, y = s4func(s3$x, deriv = 1))
s5coords <- s5func(twox[[i]][,1])

# merge y values of f"(x) segments
mergedcoords <- c(s2coords, s5coords)

# estimate constant curvature value from fitted circle
tha <- arcs[i]/radius[i] 

# plot
plot(dorsal_xsimple[[i]][,1], mergedcoords, ylim=c(-2.4, 0.1), cex=2, pch=16)
abline(h = -(tha/16), lwd=4, col='red')
title(main = paste('specimen #', i, sep=""))

# difference (nrms) in y values of the second deriv
# rms standardized by centroid size
ydis[i] <- 
  sqrt(sum((mergedcoords - (-tha/16))^2) / length(mergedcoords)) %>% 
  magrittr::divide_by(csize[i])

}

# save nrms values for comparing to circle fit in "7_compare_altmetrics.R"
write_rds(ydis, path = here('data/derived_data/compare_curvrate.rds'))
