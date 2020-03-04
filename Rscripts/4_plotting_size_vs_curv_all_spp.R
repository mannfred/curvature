library(here)
library(tidyverse)



#import data ####

# size data
size_data<-
  read.csv(
  here("data/epimedium_curv_size_data.csv"), header=TRUE) %>% 
  slice(20:77)
  

# curvature data
curv_data <-
  read.csv(
  here("data/epimedium_adj_curvature.csv"), header=TRUE) %>% 
  select(2:4)

# merge curvature data into one tibble
curv_data <-
  full_join(curv_data, size_data, ) 
  

# merge size data with curvature data ####

size_curv_data <-
  left_join(size_data, curv_data, by="species_individual_panicle_flower") %>%
  rename(species_ID = species_individual_panicle_flower) %>%
  mutate(
    species_epithet = str_extract(
      species_ID, 
      "[A-Z]+") #extacts uppercase letters from strings
  )

# plot
# For E. koreanum and E. violaceum only ####

size_curv_data<-
  left_join(size_data, curv_data, by="species_individual_panicle_flower") %>%
  rename(species_ID = species_individual_panicle_flower) %>%
  mutate(
    species_epithet = str_extract(
      species_ID, 
      "[A-Z]+") #extacts uppercase letters from strings
  ) %>% 
  filter(species_epithet == "K" | species_epithet == "V")


ggplot(
  data=size_curv_data, 
  aes(x=sepal_size_mm, y=adjusted_curvature)
) +
  geom_point(
    aes(colour=factor(species_epithet)), size=3)
