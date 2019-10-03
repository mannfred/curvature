library(here)
library(tidyverse)

#combine alignment matrix with size data ####

#alignment data
stages_days_subset <- 
  readRDS("stages_days_sort_grandiflorum.rds") %>%
  na_if("-") %>%
  drop_na() %>% #remove rows with "-"
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#size data
data5 <- readRDS("epimedium_growth_data_pivot_redefined_stages.rds")

data5_subset<- data5 %>%
  filter(Species_epithet == 'grandiflorum') %>%
  select(c(1, size)) %>%
  group_by(Species_Individual_Panicle_Flower) %>%
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  Species_Individual_Panicle_Flower, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

joiner<-
  stages_days_subset %>%
  left_join( data5_subset, by = "unique_ID")


ggplot(
  data=joiner, 
  aes(x=elapsed_days, y=size)
) +
geom_point(
  aes(x=elapsed_days, y=size), colour=group)
           
