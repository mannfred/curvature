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
# fit splines to landmarks

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


# --------------------------------
# calculate curvature from splines
# https://teazrq.github.io/stat432/rlab/spline.html
# https://stackoverflow.com/questions/35094843/get-polynomial-coefficients-from-interpolation-splines-in-r
# https://stackoverflow.com/questions/51549690/r-how-to-get-piecewise-coefficients-of-an-interpolation-spline-for-analytical-i

# test data
x1 <- seq(1, 10, 0.2)
y1 <- x1^2 
mat <- data.frame(x = x1, y =y1)

lmfit <- lm(y1 ~ splines2::naturalSpline(x1, df=4), data = mat)
plot(mat, pch=19, col='darkorange')
lines(mat$x, lmfit$fitted.values, lty=1, col='deepskyblue', lwd=4)

insMat <- naturalSpline(x1, df=4,)


# ---------
x1 <- seq(1, 10, 0.2)
y1 <- x1^2 
f <- splinefun(x1, y1, "natural") 
integrate(f, 1, 10)


f2 <- function(x) x^2
integrate(f2, 1, 10)


# ---------
x1 <- seq(1, 10, 0.2)
y1 <- x1^2 
mat <- data.frame(x = x1, y =y1)

f0 <- smooth.spline(mat)
f1predict <- predict(f0, deriv=1)

model1 <- lm(f1predict$y ~ f1predict$x)

model2func <-
  function(model) {
    
  intercept <- coef(model)[[1]]
  coeff <- coef(model)[[2]]
  
  bodyexp <- parse(text = paste(intercept, "+", coeff, "*x"))
 
  f1 <- function(x) NULL
  body(f1) <- bodyexp
  
  f2 <- Deriv::Deriv(f1)
  
  # deriv_list <- list(f1, f2)
  # return(deriv_list)
  
  kfun <- function(x) {
    f1 <- f1
    f2 <- f2
    ((f2(x))/((1 + (f1(x)^2))^1.5)) * (sqrt(1 + (f1(x))^2))
  }

  Ktot <- integrate(kfun, lower=1, upper=10)$value
  return(Ktot)
}


