library(effectsize)
library(here)
library(lmerTest)
library(Momocs)
library(tidyverse)


# --------------------------------------
# is chord length correlated with curvature in Epimedium?

# import linear measurement data
alt_metrics <- 
  read.csv(file = here("data/derived_data/alternative_metrics.csv")) %>% 
  as_tibble() 

# import total curvature data
curv_data <- 
  read_rds(here('data/derived_data/RDS_files/spline_curvature_tbl_dorsal.rds')) 

data10 <-
  data.frame(chord = alt_metrics$chord, 
             curvature = curv_data$total_K,
             arclength = curv_data$perimeter,
             species = c(rep("koreanum", 30), rep("violaceum", 27)))

ggplot(data = data10, aes(x = curvature, y = arclength, group = species, colour=species)) +
  geom_point(size=4) +
  geom_smooth(method = 'lm', se = FALSE, size=1.5) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=14)) 

# mixed effects model: t=2.016, df=54.00531, p=0.049
lmerTest::lmer(arclength ~ curvature + (1|species), data = data10) %>% summary()

#eta2 = 0.07
effectsize::t_to_eta2(-2.016, 54.00531)




# -----------------------------------------------------------------------
# is chord length correlated with curvature in for allometrically scaled curves?


# create 1000 in silico curves
coords <- vector("list", 1000)
set.seed(1)
for (i in 1:1000) {
  
  # 15 random, increasing x coordinates bw 0 and 200 
  a <- sort(sample.int(n=200, size=15))
  
  # 7 random, increasing y coordinates bw 0 and 20
  b <- sort(sample.int(n=20, size=7))
  
  # 8 random, decreasing y coordinates
  c <- c(sort(sample.int(n=20, size=8), decreasing=TRUE), b)
  
  coords[[i]] <- cbind(a,c)
  
}

# inspect curves
plot(coords[[102]])


# extract baselines
baselines_list <- 
  coords %>% 
  lapply(., function(b) c(unlist(b[,1])[1], unlist(b[,1])[15]))


# fit splines and compute curvature
curvature_tbl <-
  mapply(curvature_spline, coords, baselines_list) %>% 
  enframe() %>% 
  mutate(total_K = abs(value)*(180/pi)) 


# estimate arclength as perimeter
perim_list <- Momocs::coo_perim(Ldk(coords))


# estimate chord length
chordlength <- Momocs::coo_length(Ldk(coords))


# estimate centroid size
csize <- coo_centsize(Ldk(coords))

alldata <- 
  perim_list %>% 
  unlist() %>% 
  enframe() %>% 
  mutate(totalK = curvature_tbl$total_K,
         chord = chordlength,
         size = csize) %>% 
  rename(perimeter = value) %>% 
  mutate(acratio = perimeter/chord)

plot(alldata$perimeter, alldata$chord)
plot(alldata$perimeter, alldata$totalK)
plot(alldata$chord, alldata$totalK)

# p=0.957, t=0.054, df=998, R^2adj=-0.0009991
lm(perimeter~totalK, data=alldata) %>% summary

# p=0.836, t=-0.207, df=998, R^2adj=-0.0009589
lm(chord~totalK, data=alldata) %>% summary




# is arc:chord ratio a measure of size?
plot(alldata$acratio, alldata$size)
lm(size ~ acratio, data=alldata) %>% summary #YES. ACR is a measure of size


