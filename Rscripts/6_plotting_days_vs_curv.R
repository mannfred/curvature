library(here)
library(tidyverse)


#using days estimates from logistic modelling of size-days relationship to
#estimate days-curvature relationship

#################
#import data ####

#elapsed days/ size data

days_gran<-
  readRDS(here("data/RDS_files/size_data_grandiflorum.rds"))
  
days_kore<-
  readRDS(here("data/RDS_files/size_data_koreanum.rds"))

days_viol<-
  readRDS(here("data/RDS_files/size_data_violaceum.rds"))

#merge days/size data into one tibble

size_days_data<-
  full_join(days_gran, days_kore) %>%
  full_join(., days_viol)


###############
#curvature data
  
curv_gran<-
  read.csv(
    here("data/epimedium_curvature_grandiflorum.csv"),
    header=TRUE) %>%
  dplyr::select(4, 5) #select columns with species IDs and adjusted curvature data

curv_kore<-
  read.csv(
    here("data/epimedium_curvature_koreanum.csv"),
    header=TRUE) %>%
  dplyr::select(4, 5)

curv_viol<-
  read.csv(
    here("data/epimedium_curvature_violaceum.csv"),
    header=TRUE) %>%
  dplyr::select(4, 5)


#####################################
#merge curvature data into one tibble

curv_data<-
  full_join(curv_gran, curv_kore) %>%
  full_join(., curv_viol)


#merge curvature data with size/days data

all_data<-
  left_join(size_days_data, curv_data, by="species_individual_panicle_flower") %>%
  mutate(
     species_epithet = str_extract(
       species_individual_panicle_flower, 
      "[A-Z]+") #extacts uppercase letters from strings
  )

#####
#plot

ggplot(
  data=all_data, 
  aes(x=predicted_days_elapsed, y=adjusted_curvature)
) +
  geom_point(
    aes(colour=factor(species_epithet)), size=2.1)
