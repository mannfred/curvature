library(here)
library(tidyverse)

#for E. grandiflorum only (all 3 spp. done simultaneously below)####

#combine alignment matrix with size data 

#alignment data
stages_days_subset <- #only import E. grandiflorum data
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

joiner_grandiflorum<-
  stages_days_subset %>%
  left_join( data5_subset, by = "unique_ID")

#plot single species
ggplot(
  data=joiner_grandiflorum, 
  aes(x=elapsed_days, y=size)
) +
  geom_point(
    aes(x=elapsed_days, y=size))


##########################################
# make scatterplot with all 3 species ####

#import elapsed days data
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

#import size data
data5 <- readRDS("epimedium_growth_data_pivot_redefined_stages.rds")

size_data <- data5 %>%
  select(c(1, size)) %>% #choose columns with species IDs and size data
  group_by(Species_Individual_Panicle_Flower) %>%
  mutate(row_ID = row_number()) %>% #create unique identifiers by row number
  unite(unique_ID,  Species_Individual_Panicle_Flower, row_ID, sep="_") %>% #create unique_ID column
  ungroup()


#for all three spp.
joiner<-list(size_data, gran_stages_days, kore_stages_days, viol_stages_days) %>%
  reduce(full_join, by = "unique_ID") %>%
  mutate(days_elapsed = coalesce(elapsed_days.x, elapsed_days.y, elapsed_days)) %>% #create merged days column
  mutate(flower_stage = coalesce(stage.x, stage.y, stage)) %>% #create merged stage column
  select(-c(3:8)) %>% #remove old columns
  mutate(
    species_epithet = str_extract(
      unique_ID, 
      "[a-z]+") #extacts lowercase letters from strings
  ) %>%
  group_by(species_epithet)

#plot
ggplot(
  data=joiner, 
  aes(x=days_elapsed, y=size)
) +
  geom_point(
    aes(colour=factor(species_epithet)))

