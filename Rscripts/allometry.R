library(here)
library(Momocs)
library(pracma)
library(tidyverse)


# allometry



# --------------------------------------------------------------
# is chord length and arc length correlated for the unit circle?

# a set of five segments of the unit circle
# c1 = sqrt(3)/2
# c2 = 1,0 to sqrt(2)/2, sqrt(2)/2
# c3 = 1,0 to 0,1
# 

f1 <- function(t) c(cos(t), sin(t))
f3 <- function(theta) 2*sin(theta/2)
x1 <- seq(0, 1, by=0.01)
y1 <- f3(x1)

plot(x1,y1)

# for circles, arc:chord and mand. index measure curvature 
# independent of size

f1 <- function(t) c(cos(t), sin(t))
f2 <- function(t) c(2*cos(t), 2*sin(t))
f3 <- function(t) c(3*cos(t), 3*sin(t))
arclength(f1, 0, pi) #3.14
arclength(f2, 0, pi) #6.28
arclength(f3, 0, pi) #9.42





# ----------------------------------------------
# for polynomials
# Large specimen: arcchord ratio = 1.027433
f3 <- function(t) c(t, t^2)
arclength(f3, 0, 3) #9.747089
ed(c(0,0), c(3,9)) #9.486833


# Small specimen: arcchord ratio = 1.027433
f4 <- function(t) c(t, 2*t^2)
arclength(f4, 0, 1.5) #4.873544
ed(c(0,0), c(1.5,4.5)) #4.743416


# tiny specimen: arcchord ratio = 1.027433
f5 <- function(t) c(t, 4*t^2)
arclength(f5, 0, 0.75) #2.436772
ed(c(0,0), c(0.75,2.25)) #2.371708




