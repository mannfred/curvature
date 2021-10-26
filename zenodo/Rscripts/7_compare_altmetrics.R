# compare historical metrics to total and point-wise curvature 


library(effectsize)
library(emmeans)
library(GGally)
library(Hmisc)
library(lmerTest)
library(tidyverse)
library(here)


# ------------------------------------------
# import alt curv metrics measured in "6_measure_altmetrics.R"
# import total curvature measured in "3_curvature_splines_dorsal.R"
# import circles fitted in "3_circle_fitting.R" 
# import point-curvature measured in "3_pointwise_k.R" 


# import historic metrics ("6_measure_altmetrics.R")
alt_metrics <- 
  read.csv(file = here("data/derived_data/alternative_metrics.csv")) %>% 
  as_tibble() 

# import total curvature data ("3_curvature_splines_dorsal.R")
curv_data <- 
  read_rds(here('data/derived_data/RDS_files/spline_curvature_tbl_dorsal.rds')) 

# import pracma circles (from "3_circle_fitting.R")
circledata <- read.csv(file=here('data/derived_data/circles_fitby_pracma.csv'))

# import curvature rate comparisons ("3_pointwise_k.R")
ydis <- read_rds(path=here('data/derived_data/compare_curvrate.rds'))



# ---------------------------------------
# define variables

totalK <- curv_data$total_K
angle_of_deflection <- alt_metrics$theta
arc_chord_ratio <- circledata$arcs / alt_metrics$chord
chord_versine_ratio <- alt_metrics$chord / alt_metrics$depth
iradius <- (1/circledata$radius)*circledata$arcs #to adjust 1/R by arclength
csize <- circledata$csize
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
         species,
         ydis) %>% 
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



# pairwise comparisons of metrics for E. koreanum
pairK <- Hmisc::rcorr(as.matrix(all_metrics)[1:30, 1:5], type="spearman") 


# write.csv(data.frame(pairK$r), file=here('data/pairwiseK_r.csv'))
# write.csv(data.frame(pairK$P), file=here('data/pairwiseK_p.csv'))


# pairwise comparisons of metrics for E. violaceum
pairV <- Hmisc::rcorr(as.matrix(all_metrics)[31:57, 1:5], type="spearman")

# write.csv(data.frame(pairV$r), file=here('data/pairwiseV_r.csv'))
# write.csv(data.frame(pairV$P), file=here('data/pairwiseV_p.csv'))


# rho values associated with totalK vs not associated with totalK
tk <- as.numeric(c(pairK$r[2:5,1], pairV$r[2:5,1]))
ot <- as.numeric(c(pairK$r[3:5,2], pairK$r[4:5,3], pairK$r[5,4], pairV$r[3:5,2], pairV$r[4:5,3], pairV$r[5,4]))
fc <- as.factor(c(replicate(8, "totalk"), replicate(12, "other")))
tkot <- data.frame(rho = abs(c(tk, ot)), method = fc)

lm(rho ~ method, data=tkot) %>% 
  emmeans(., specs = "method") %>% 
  pairs()


# -----------------------------------------
# import circlefit estimates: this data was `print()`ed by pracma::circlefit() and pasted into this .csv

# lower RMS error is a better circle fit
rms_error <- 
  read.csv(file=here('data/derived_data/circle_rms_error_pracma.csv')) %>% 
  rename(error = ï..rms_error) %>%
  mutate(norm_error = error/circledata$csize) # normalize rms error by centroid size
  


# ---------------------------------
# test for correlation bw circle fit  
# and nrms of point-curvature

comparedata <- 
  data.frame(
    norm_error = rms_error$norm_error, 
    ydis = ydis,
    species = species) 


mytheme <- 
  theme_bw() + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))


ggplot(data = comparedata, aes(x = norm_error, y = ydis, colour=species)) +
  geom_point(aes(size = 2.5, alpha=1.5)) +
  geom_smooth(method = 'lm', se = FALSE, size=1.5) +
  scale_colour_manual(values=unique(colour_ids)) +
  # geom_text(aes(label=curv_data$name)) +
  mytheme


# is circle fit correlated with disparity bw Ktot and 1/R metrics?
lmerTest::lmer(ydis ~ norm_error + (1|species), data = comparedata) %>% summary()

# effect size (partial eta2) = 0.19
# p = 0.000728
effectsize::t_to_eta2(3.580, 54.88682)


