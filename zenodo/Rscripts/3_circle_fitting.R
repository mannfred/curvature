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


# -------------------------
# fit circles to landmarks

# estimate radius of circle
# extract columns 1 and 2 (x and y) from a list of matrices
xcoords <- lapply(dorsal_x$coo, function(i) i[,1])
ycoords <- lapply(dorsal_x$coo, function(i) i[,2])

# fit circles 
# RMS error is printed to console, 
# so it has been saved to "circle_rms_error_pracma.csv"
circles <- mapply(pracma::circlefit, xcoords, ycoords)

# extract radii
radius <- circles[3,]


# bounded arclength for fitted circles to adjust 1/R metric
# see: "7_compare_altmetrics.R"
arcs <- vector(); r0 <- numeric(); f <- list()
a1 <- numeric(); a2 <- numeric(); 
v1 <- list(); v2 <- list(); angle <- numeric()

for (i in 1:length(dorsal_x$coo)) {
  # radii
  r0[i] <- as.numeric(circles[3,i])
  
  # parametrized circle functions
  f[[i]] <- function(t) c(r0[i] * cos(t), r0[i] * sin(t))
  
  # the subtraction of two points produces a displacement vector between the two points
  #  https://math.stackexchange.com/questions/873631/subtracting-two-points-the-math
  v1[[i]] <- dorsal_x$coo[[i]][1,] -  circles[1:2,i]  #vector from circle's center to first x
  v2[[i]] <- dorsal_x$coo[[i]][16,] - circles[1:2,i]  #vector from circle's center to last x
  
  # angle between xaxis and v1 and v2
  a1[i] <- atan2(v1[[i]][2], v1[[i]][1])
  a2[i] <- atan2(v2[[i]][2], v2[[i]][1])
  
  # angle bw v1 and v2
  angle[i] <- a1[i] - a2[i]
  
  # feed angle start-end positions to arclength function
  arcs[i] <- arclength(f[[i]], 0, angle[i])$length 
}

# compute centroid size to normalize rms error estimates
# see: "7_compare_altmetrics.R"
csize <- coo_centsize(dorsal_x)

# compile circledata 
circledata <-  cbind(radius, arcs, csize)

# save to .csv for use in "7_compare_altmetrics.R"
# write.csv(circledata, file=here('data/derived_data/circles_fitby_pracma.csv'))



# -------------------------------------
# plotting circles
for (i in 1:length(circles)){
  
  plot(dorsal_xsimple[[i]])
  title(main = paste('specimen #', i, sep=""))
  
  x0 <- circles[1,i]
  y0 <- circles[2,i]
  r0 <- circles[3,i]
  
  w  <- seq(0, 2*pi, len=100)
  xx <- r0 * cos(w) + x0
  yy <- r0 * sin(w) + y0
  
  lines(xx, yy, col="red", lw=2)
}
