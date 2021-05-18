library(here)
library(Momocs)
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

# linear measurements
depth <- coo_width(dorsal_x) %>%  as.numeric()
chord <- coo_length(dorsal_x) %>% as.numeric()

# angle of declension
# angle bw two vectors: https://stackoverflow.com/questions/21483999/using-atan2-to-find-angle-between-two-vectors
angles <- numeric()
for (i in 1:length(dorsal_x)) {

center <- coo_centpos(dorsal_x$coo[[i]])

a <- center - dorsal_x$coo[[i]][1,] #vector from apex of spur to centroid
b <- center - dorsal_x$coo[[i]][15,] #vector from base of spur to centroid
angles[i] <- 180 - (180/pi) * (atan2(b[2], b[1]) - atan2(a[2], a[1]))
}


# compile and write to .csv
alt_metrics <- data.frame(chord = chord, depth = depth, theta = angles)
# write.csv(alt_metrics, file=here('data/derived_data/alternative_metrics.csv'))
