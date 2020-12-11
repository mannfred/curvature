library(tidyverse)
library(curvr)
library(Momocs)

# construct data
x <- c(1, sqrt(3)/2, sqrt(2)/2, 1/2, 0)
y <- c(0, 1/2, sqrt(2)/2, sqrt(3)/2, 1)

cdat <-
  cbind(x, y) %>%
  list() %>%
  Ldk() %>%
  coo_alignxax()

# fit polynomial
cpoly <- npoly(cdat$coo, degree = 4)

# plot LMs
coo_plot(cdat$coo)

# draw lines
cpoly$shp1 %>%
  npoly_i() %>%
  coo_draw(border='red')


#compute Ktot
total_curvature(cpoly$shp1, c(-0.2588190, 0.2588190))

