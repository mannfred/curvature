library(effectsize)
library(emmeans)
library(GGally)
library(tidyverse)
library(here)


# ------------------------------------------
# import alt curv metrics measured in "6_measure_altmetrics.R"

alt_metrics <- 
  read.csv(file = here("data/derived_data/alternative_metrics.csv")) %>% 
  as_tibble() 

# import total curvature data
curv_data <- 
  read_rds(here('data/derived_data/RDS_files/spline_curvature_tbl_dorsal.rds')) 

# import pracma circles (from "3_curvature_splines_dorsal.R")
pracma_circles <- read.csv(file=here('data/derived_data/circles_fitby_pracma.csv'))



#create vars
totalK <- curv_data$total_K
angle_of_deflection <- alt_metrics$theta
arc_chord_ratio <- curv_data$perimeter / alt_metrics$chord
chord_versine_ratio <- alt_metrics$chord / alt_metrics$depth
# radius <- ((alt_metrics$chord / 2) / sin(angle_of_deflection*(pi/180))) #sin function operates on radians
iradius <- (1/pracma_circles[,2])*curv_data$perimeter #to adjust 1/R by arclength
species <- as.factor(c(replicate(30, "koreanum"), replicate(27, "violaceum")))

# -------------------------
# plot pairwise comparisons

#colours from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 
colour_ids <- c("#009E73", "#CC79A7")

group_ids <- c("koreanum", "violaceum")

all_metrics <-
  tibble(totalK, 
         angle_of_deflection, 
         arc_chord_ratio, 
         chord_versine_ratio, 
         iradius,
         species) %>% 
  group_by(species)

#store plot
pairwise_plot <-
  ggpairs(data = all_metrics, columns = 1:5, mapping = aes(colour = species),  
        columnLabels = c("total curvature", "angle of deflection", "arc:chord ratio", "mandibular index", "inverse radius"),
        
        lower = list(continuous = function(data, mapping, ...) {
          ggally_smooth_loess(data = data, mapping = mapping, alpha =0.5, size=3.5, se=F) + 
            scale_colour_manual(values=colour_ids)}),
        
        diag  = list(continuous = function(data, mapping, ...) {
          ggally_densityDiag(data = data, mapping = mapping, alpha = 0.7) +
            scale_fill_manual(values=colour_ids, limits=group_ids)}),
        
        upper = list(continuous = 'blank')
)

# show plot
pairwise_plot + theme_classic()
 



# ----------------------------------
# run stats on pairwise comparisons

library(Hmisc)

# pairwise comparisons of metrics for E. koreanum
pairK <- Hmisc::rcorr(as.matrix(all_metrics)[1:30, 2:6], type="spearman") 


# write.csv(data.frame(pairK$r), file=here('data/pairwiseK_r.csv'))
# write.csv(data.frame(pairK$P), file=here('data/pairwiseK_p.csv'))


# pairwise comparisons of metrics for E. violaceum
pairV <- Hmisc::rcorr(as.matrix(all_metrics)[31:57, 2:6], type="spearman")

# write.csv(data.frame(pairV$r), file=here('data/pairwiseV_r.csv'))
# write.csv(data.frame(pairV$P), file=here('data/pairwiseV_p.csv'))



# ---------------------------------
# inspect outliers

plot(totalK, iradius, text(x=totalK, iradius, labels=curv_data$name))

# koreanum outliers:
# 1_6_2(4), 1_8_3(10),  2_3_2(13), 2_5_1(16),
# 2_9_2(21), 2_9_3(22), 2_10_1(23)

# violaceum outliers:
# 2_3_5(44), 2_3_3(42)

# residuals between total_curvature and 1/R
res <- lm(iradius ~ totalK)$residuals 
plot(1:57, res, text(x=1:57, res, labels=curv_data$name))


# this data `print()`ed by pracma::circlefit() and pasted into this .csv
# lower RMS error is a better circle fit
rms_error <- 
  read.csv(file=here('data/derived_data/circle_rms_error_pracma.csv')) %>% 
  rename(rms_error = ï..rms_error)
  
plot(1:57, rms_error$rms_error, text(x=1:57, rms_error$rms_error, labels=curv_data$name))


# ---------------------------------
# plot 

circledata <- 
  data.frame(
    rms_error = rms_error$rms_error, 
    residuals = abs(res), 
    species = c(rep("koreanum", 30), rep("violaceum", 27))) %>% 
  group_by(species)


mytheme <- 
  theme_bw() + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))


ggplot(data = circledata, aes(x = rms_error, y = residuals, group=species, color=species)) +
  geom_point(aes(size = 2.5, alpha=1.5)) +
  geom_smooth(method = 'lm', se = FALSE, size=1.5) +
  scale_colour_manual(values=unique(colour_ids)) +
  # geom_text(aes(label=curv_data$name)) +
  mytheme

# do shapes fitted by circles agree more with Ktot?
lm(residuals ~ rms_error*species, data=circledata) %>% summary()
# t = 7.003, df_error = 53

# effect size (partial eta2) = 0.48
effectsize::t_to_eta2(7.003, 53)



