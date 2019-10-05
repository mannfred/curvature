library(here)
library(tidyverse)

#logistic models
#for E. grandiflorum ####

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

# logistic model parameters estimated by max lik ####

# extract max value to normalize size range to (0,1)
joiner_grandiflorum$adjusted_size<-
  joiner_grandiflorum$size/max(joiner_grandiflorum$size) 


model <- glm(adjusted_size ~ elapsed_days, joiner_grandiflorum, family = quasibinomial)

#extrapolate from model
#generate x-range
model_df <- data.frame(elapsed_days = seq(0, 30, 1))

#predict the size values (as a probability) using the above model
model_df$size <- predict(model, newdata=model_df, type="response")

# plot the modeled probability values
ggplot(model_df, aes(x=elapsed_days, y=size)) + geom_line()




################################################
#plotting data points with regression curve ####
ggplot(joiner_grandiflorum, aes(x=elapsed_days, y=adjusted_size)) + 
  geom_point(alpha=.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, method.args = list(family=quasibinomial)) + 
  ylab("size (mm)") + xlim(0,21)

#assesing model fit with pseudo R2
# > model$null.deviance
# [1] 159.2194
# > model$deviance
# [1] 40.83834

model_chi <- model$null.deviance - model$deviance
pseudo_r2 <- model_chi / model$null.deviance
#pseudo R2 = 0.7435

#assesing model significance (approximating p value)
chi_df <- model$df.null - model$df.residual

chisq_prob <- 1 - pchisq(model_chi, chi_df)
#chisq_prob = 0 



############################################################
# using model to predict elapsed days from curvature dataset

#get the equation from this model to plug in "size" and get "elapsed_days"

size_data<-
  read.csv(
    here("data/epimedium_curv_size_data.csv"), 
    header=TRUE) %>%
  slice(1:19) #select grandiflorum only

size_data$adjusted_size<-
  size_data$sepal_size_mm/max(size_data$sepal_size_mm) 


size_data$elapsed_days<- 
  predict(model, newdata=size_data$adjusted_size, type="response")
