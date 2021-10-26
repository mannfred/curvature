# measure curvature as defined by methods 2-5 (see: main text)
# need to measure several linear features of each specimen: chord, versine, angle of deflection

library(geomorph)
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
dorsal <- epi_lmk_data[c(31, 1:15),,]  %>%  Momocs::Ldk()

# align the longest axis of each shape along the x-axis
dorsal_x <- Momocs::coo_alignxax(dorsal)

# linear measurements
depth <- coo_width(dorsal_x) %>%  as.numeric()
chord <- coo_length(dorsal_x) %>% as.numeric()

# angle of declension
# angle bw two vectors: https://stackoverflow.com/questions/21483999/using-atan2-to-find-angle-between-two-vectors
angles <- numeric(); angles2 <- numeric() 
v1 <- list(); v2 <- list()
a1 <- numeric(); a2 <- numeric()

for (i in 1:length(dorsal_x)) {

v1[[i]] <- dorsal_x$coo[[i]][1,] - dorsal_x$coo[[i]][8,]  #vector from base of spur to midpoint
v2[[i]] <- dorsal_x$coo[[i]][1,] - dorsal_x$coo[[i]][16,] #vector from base of spur to apex

# angle between xaxis and v1 and v2
a1[i] <- atan2(v1[[i]][2], v1[[i]][1])
a2[i] <- atan2(v2[[i]][2], v2[[i]][1])

angles[i] <- (a1[i]-a2[i])*(180/pi)

angles2[i] <- 
  
  if(angles[i] < 0){ 
  angles[i] + 360
    
  } else {angles[i]}

}



# compile and write to .csv
alt_metrics <- data.frame(chord = chord, depth = depth, theta = angles2)
write.csv(alt_metrics, file=here('data/derived_data/alternative_metrics.csv'))
