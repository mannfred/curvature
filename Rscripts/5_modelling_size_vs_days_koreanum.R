library(here)
library(tidyverse)
library(drc) #for logistic models 

#logistic models
#for E. koreanum ####

#######################################################
#combine alignment matrix (elapsed days) with size data 

#import alignment data
stages_days_kore <- #import E. koreanum data
  readRDS(here("data/RDS_files/stages_days_sort_koreanum.rds")) %>%
  na_if("-") %>% 
  drop_na() %>% #remove rows with "-" 
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#import size data
data5 <- readRDS(here("data/RDS_files/epimedium_growth_data_pivot_redefined_stages.rds"))

data5_kore <- data5 %>%
  filter(Species_epithet == 'koreanum') %>%
  dplyr::select(c(1, size)) %>% #select masked by drc/aomisc packages 
  group_by(Species_Individual_Panicle_Flower) %>%
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  Species_Individual_Panicle_Flower, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#join alignment and size data
joiner_koreanum<-
  stages_days_kore %>%
  left_join( data5_kore, by = "unique_ID")


#plotting elapsed_days vs size 
ggplot(joiner_koreanum, aes(x=elapsed_days, y=size)) + 
  geom_point(alpha=.5) +
  ylab("size (mm)") + xlim(0,27)


#################################################################
# using model to predict elapsed days from curvature dataset ####

model_kore <-
  drm(size ~ elapsed_days, fct = L.4(), data=joiner_koreanum)

plot(model_kore, log="", main = "Logistic function")

#import size data used to make curvature measurements
size_data_kore<-
  read.csv(
    here("data/epimedium_curv_size_data.csv"), 
    header=TRUE) %>%
  slice(20:50) #select koreanum only

#use model to predict "elapsed days" from size data
size_data_kore$predicted_days_elapsed<-
  predict(model_kore, data.frame(size_data_kore$sepal_size_mm))

saveRDS(size_data_kore, file=here("data/RDS_files/size_data_koreanum.rds"))

#plot 
ggplot(size_data_kore, aes(x=sepal_size_mm, y=predicted_days_elapsed)) + 
  geom_point(alpha=.5) +
  ylab("predicted_elapsed_days") + xlim(0,35)


###########################
#append curvature data ####

#import curvature data
curv_kore<-
  read.csv(
    here("data/epimedium_curvature_koreanum.csv"),
    header=TRUE) %>%
  dplyr::select(4, 5)


