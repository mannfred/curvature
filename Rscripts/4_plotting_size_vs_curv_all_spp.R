library(here)
library(tidyverse)



#import data ####

# size data
size_data<-
  read.csv(
  here("data/epimedium_curv_size_data.csv"), header=TRUE) %>% 
  slice(20:77) %>% 
  rename(name = species_individual_panicle_flower)
  

# curvature data
curv_data <-
  read.csv(
  here("data/epimedium_adj_curvature.csv"), header=TRUE) %>% 
  select(2:5)

# merge curvature data into one tibble
curv_size_data <-
  left_join(curv_data, size_data, by ='name') %>% 
  as_tibble() %>% 
  mutate(species = str_extract(name, "[A-Z]+")) # extacts uppercase letters from strings
  
  
# plot sepal size vs curv
ggplot(
  data=curv_size_data, 
  aes(x=sepal_size_mm, y=adjusted_curvature)) +
geom_point(
  aes(colour=factor(species)), size=3) +
labs(x = "Sepal length (mm)", 
     y = "Total mean curvature (degrees/mm)",
     colour = "species") +
scale_colour_manual(name= "Species", 
                    breaks = c("K", "V"), 
                    labels = c(expression(italic("E. koreanum")), 
                               expression(italic("E. violaceum"))),
                    values = c("#009E73", "#CC79A7")) +
theme_bw() + 
theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      legend.position=c(.85, .65))
               

# plot sepal size vs arc length
ggplot(
  data=curv_size_data, 
  aes(x=sepal_size_mm, y=arclength)) +
  geom_point(
    aes(colour=factor(species)), size=3) +
  labs(x = "Sepal length (mm)", 
       y = "Arc length (mm)",
       colour = "species") +
  scale_colour_manual(name= "Species", 
                      breaks = c("K", "V"), 
                      labels = c(expression(italic("E. koreanum")), 
                                 expression(italic("E. violaceum"))),
                      values = c("#009E73", "#CC79A7")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position=c(.85, .25))

