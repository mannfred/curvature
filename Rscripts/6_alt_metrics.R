library(tidyverse)
library(here)


# ------------------------
# import alt curv metrics
alt_metrics <- 
  read.csv(file = here("data/alternative_metrics_comparison.csv")) %>% 
  as_tibble() %>% 
  mutate(theta = 360 - theta)


# import total curvature data
curv_data <- 
  read.csv(here('data/epimedium_adj_curvature_dorsal.csv')) 

# import pracma circles
pracma_circles <- read.csv(file=here('data/pracma_circles.csv'))

#create vars
totalK <- curv_data$total_curvature
angle_of_deflection <- alt_metrics$theta
arc_chord_ratio <- curv_data$arclength / alt_metrics$chord
chord_versine_ratio <- alt_metrics$chord / alt_metrics$depth
# radius <- ((alt_metrics$chord / 2) / sin(angle_of_deflection*(pi/180))) #sin function operates on radians
radius <- pracma_circles[,2]


# -------------------------
# plot pairwise comparisons

#colours from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 
colour_ids <-
  c(replicate(30, "#E69F00"), #31 E. koreanum samples
    replicate(27, "#56B4E9")) #27 E. violaceum samples


#pairwise plots
library(GGally)

all_metrics <-
  tibble(colour_ids,
         totalK, 
         angle_of_deflection, 
         arc_chord_ratio, 
         chord_versine_ratio, 
         1/radius)

#store plot
pairwise_plot <- 
  ggpairs(
    all_metrics,
    columns = 2:6,
    mapping = aes(colour = colour_ids, alpha=0.35),
    columnLabels = 
      c("total curvature",
        "angle of deflection",
        "arc:chord ratio",
        "mandibular index",
        "inverse radius"),
    lower = list(continuous = 'smooth', alpha = 0.5),
    diag  = list(continuous = 'densityDiag', alpha = 0.3),
    upper = list(continuous = 'blank'))
    

#show plot
pairwise_plot + theme_classic()



# ----------------------------------
# run stats on pairwise comparisons

library(Hmisc)

# pairwise comparisons of metrics for E. koreanum
pairK <- Hmisc::rcorr(as.matrix(all_metrics)[1:30, 2:6], type="pearson") 


# write.csv(data.frame(pairK$r), file=here('data/pairwiseK_r.csv'))
# write.csv(data.frame(pairK$P), file=here('data/pairwiseK_p.csv'))


# pairwise comparisons of metrics for E. violaceum
pairV <- Hmisc::rcorr(as.matrix(all_metrics)[31:57, 2:6], type="pearson")

# write.csv(data.frame(pairV$r), file=here('data/pairwiseV_r.csv'))
# write.csv(data.frame(pairV$P), file=here('data/pairwiseV_p.csv'))

# ----------------------
# inspect outliers

plot(1/radius, 
     arc_chord_ratio, 
     text(x=1/radius, 
          y=arc_chord_ratio, 
          labels=curv_data$name))

# outliers:
# 1_8_3(10), 2_3_2(13), 2_3_3(14), 2_9_2(21),  2_9_3(22), 2_10_1(23) 


# new outliers 
# rms K233 (14) = 0.1792714
# rms K292(21) = 0.4432476
# rms V235(44) = 0.6492074 

# lower RMS error is a better circle fit
rms_error <- 
  read.csv(file=here('data/circle_rms_error.csv')) 
  
# root of squared residuals between total_curvature and 1/R
res <- sqrt(residuals(lm(totalK~ 1/radius))^2)

# plot 
plot(rms_error$rms_error, 
     res,
     text(x=rms_error$rms_error, 
          y=res, 
          labels=curv_data$name))

summary(lm(rms_error$rms_error ~ res))
 