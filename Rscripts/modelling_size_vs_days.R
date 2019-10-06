library(here)
library(tidyverse)
library(drc) #for logistic models 

#logistic models
#for E. grandiflorum ####

#combine alignment matrix (elapsed days) with size data 

#import alignment data
stages_days_subset <- #only import E. grandiflorum data
  readRDS("stages_days_sort_grandiflorum.rds") %>%
  na_if("-") %>% 
  drop_na() %>% #remove rows with "-" 
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#import size data
data5 <- readRDS("epimedium_growth_data_pivot_redefined_stages.rds")

data5_subset<- data5 %>%
  filter(Species_epithet == 'grandiflorum') %>%
  dplyr::select(c(1, size)) %>% #select masked by drc/aomisc packages 
  group_by(Species_Individual_Panicle_Flower) %>%
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  Species_Individual_Panicle_Flower, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#join alignment and size data
joiner_grandiflorum<-
  stages_days_subset %>%
  left_join( data5_subset, by = "unique_ID")


#plotting elapsed_days vs size 
ggplot(joiner_grandiflorum, aes(x=elapsed_days, y=adjusted_size)) + 
  geom_point(alpha=.5) +
  ylab("size (mm)") + xlim(0,21)




#################################################################
# using model to predict elapsed days from curvature dataset ####

model <-
  drm(size ~ elapsed_days, fct = L.4(), data=joiner_grandiflorum)

plot(model, log="", main = "Logistic function")

#import size data used to make curvature measurements
size_data<-
  read.csv(
    here("data/epimedium_curv_size_data.csv"), 
    header=TRUE) %>%
  slice(1:19) #select grandiflorum only

#use model to predict "elapsed days" from size data
size_data$predicted_days_elapsed<-
  predict(model, data.frame(size_data$sepal_size_mm))

ggplot(size_data, aes(x=sepal_size_mm, y=predicted_days_elapsed)) + 
  geom_point(alpha=.5) +
  ylab("predicted_elapsed_days") + xlim(0,35)




