library(here)
library(tidyverse)

#import data ####

#size data
size_data<-
  read.csv(
  here("data/epimedium_curv_size_data.csv"), 
  header=TRUE)

#curvature data

curv_gran<-
  read.csv(
    here("data/epimedium_curvature_grandiflorum.csv"),
    header=TRUE) %>%
  select(4, 5) #select columns with species IDs and adjusted curvature data
  
curv_kore<-
  read.csv(
    here("data/epimedium_curvature_koreanum.csv"),
    header=TRUE) %>%
  select(4, 5)

curv_viol<-
  read.csv(
    here("data/epimedium_curvature_violaceum.csv"),
    header=TRUE) %>%
  select(4, 5)

#merge curvature data into one tibble

curv_data<-
  full_join(curv_gran, curv_kore) %>%
  full_join(., curv_viol)

# merge size data with curvature data ####

size_curv_data<-
  left_join(size_data, curv_data, by="species_individual_panicle_flower") %>%
  rename(species_ID = species_individual_panicle_flower) %>%
  mutate(
    species_epithet = str_extract(
      species_ID, 
      "[A-Z]+") #extacts uppercase letters from strings
  )

#plot

ggplot(
  data=size_curv_data, 
  aes(x=sepal_size_mm, y=adjusted_curvature)
) +
  geom_point(
    aes(colour=factor(species_epithet)), size=3)
