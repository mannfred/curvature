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
  aes(x=days_elapsed, y=size)
) +
geom_point(
  aes(colour=factor(species_epithet)))

########### make scatterplot w 3 species

gran_stages_days <- 
  readRDS("stages_days_sort_grandiflorum.rds") %>%
  na_if("-") %>%
  drop_na() %>% #remove rows with "-"
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

kore_stages_days <-
  readRDS("stages_days_sort_koreanum.rds") %>%
  na_if("-") %>%
  drop_na() %>% #remove rows with "-"
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

viol_stages_days <-
  readRDS("stages_days_sort_violaceum.rds") %>%
  na_if("-") %>%
  drop_na() %>% #remove rows with "-"
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#size data
data5 <- readRDS("epimedium_growth_data_pivot_redefined_stages.rds")

data5_subset<- data5 %>%
  # filter(Species_epithet == 'grandiflorum') %>%
  select(c(1, size)) %>%
  group_by(Species_Individual_Panicle_Flower) %>%
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  Species_Individual_Panicle_Flower, row_ID, sep="_") %>% #create unique_ID column
  ungroup()



joiner<-list(data5_subset, gran_stages_days, kore_stages_days, viol_stages_days) %>%
  reduce(full_join, by = "unique_ID") %>%
  mutate(days_elapsed = coalesce(elapsed_days.x, elapsed_days.y, elapsed_days)) %>% #create merged days column
  mutate(flower_stage = coalesce(stage.x, stage.y, stage)) %>% #create merged stage column
  select(-c(3:8)) %>% #remove old columns
  mutate(
    species_epithet = str_extract(
      unique_ID, 
      "[a-z]+") #extacts words from strings
  ) %>%
  group_by(species_epithet)


