library(here)
library(tidyverse)
library(drc) #for logistic models 

#logistic models
#for E. violaceum ####

#######################################################
#combine alignment matrix (elapsed days) with size data 

#import alignment data
stages_days_viol <- #import E. violaceum data
  readRDS(here("data/RDS_files/stages_days_sort_violaceum.rds")) %>%
  na_if("-") %>% 
  drop_na() %>% #remove rows with "-" 
  group_by(ID) %>% #create groups by flower
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  ID, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#import size data
data5 <- readRDS(here("data/RDS_files/epimedium_growth_data_pivot_redefined_stages.rds"))

data5_viol <- data5 %>%
  filter(Species_epithet == 'violaceum') %>%
  dplyr::select(c(1, size)) %>% #select masked by drc/aomisc packages 
  group_by(Species_Individual_Panicle_Flower) %>%
  mutate(row_ID = row_number()) %>%
  unite(unique_ID,  Species_Individual_Panicle_Flower, row_ID, sep="_") %>% #create unique_ID column
  ungroup()

#join alignment and size data
joiner_violaceum<-
  stages_days_viol %>%
  left_join( data5_viol, by = "unique_ID")


#plotting elapsed_days vs size 
ggplot(joiner_violaceum, aes(x=elapsed_days, y=size)) + 
  geom_point(alpha=.5) +
  ylab("size (mm)") + xlim(0,22)



#################################################################
# using model to predict elapsed days from curvature dataset ####

model_viol <-
  drm(size ~ elapsed_days, fct = L.4(), data=joiner_violaceum)

plot(model_viol, log="", main = "Logistic function")

#import size data used to make curvature measurements
size_data_viol<-
  read.csv(
    here("data/epimedium_curv_size_data.csv"), 
    header=TRUE) %>%
  slice(51:77) #select violaceum only

#use model to predict "elapsed days" from size data
size_data_viol$predicted_days_elapsed<-
  predict(model_viol, data.frame(size_data_viol$sepal_size_mm)) #predict function from package 'drc'

saveRDS(size_data_viol, file=here("data/RDS_files/size_data_violaceum.rds"))

#plot 
ggplot(size_data_viol, aes(x=sepal_size_mm, y=predicted_days_elapsed)) + 
  geom_point(alpha=.5) +
  ylab("predicted_elapsed_days") + xlim(0,35)

