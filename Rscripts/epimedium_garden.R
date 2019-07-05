######################

#power analysis for epimedium exp.

pwr.t.test(d=0.30, sig.level = 0.05, power=0.7, type="two.sample", alternative = 'greater')

#Two-sample t test power calculation 

#n = 105.253 
#d = 0.3
#sig.level = 0.05
#power = 0.7
#alternative = greater

#NOTE: n is number in *each* group
#........need 105 samples of each epimedium spp.......



#plotting garden data####

library(here)
library(stringr) 
library(tidyverse) #dplyr, ggplot2


#pivot original data frame so that "stage" info can be manually inputted 

data<-
  read.csv('epimedium_growth_data.csv', header=TRUE) %>%
  select(Species_Individual_Panicle_Flower, 2:25) %>% #add ", 2:25" parameter to select data upto column 25 (col25=May 2, 2019)
  as.tibble() %>%
  gather(key="date", value="size", -Species_Individual_Panicle_Flower)

levels<- unique(data$Species_Individual_Panicle_Flower) %>%
  str_sort(., numeric = TRUE) #create a levels vector with the identifiers in the right order (stringr)

data_sort<-
  data %>%
  mutate(Species_Individual_Panicle_Flower = factor(Species_Individual_Panicle_Flower, levels=levels)) %>%
  arrange(Species_Individual_Panicle_Flower) %>% 
  print(n=50) %>%
  write.csv("epimedium_growth_data_pivot.csv", row.names=F) #stage info will be manually entered to this file :(


#visualize stage-size relationship####

data2<-
  read.csv('epimedium_growth_data_pivot.csv', header=TRUE) %>%
  filter(grepl('koreanum', Species_Individual_Panicle_Flower)) %>% #change species name to generate separate graphs
  na.omit() %>%
  mutate(stage = factor(stage, levels=c("E", "G", "O", "P", "D", "A"))) #reorder $stage 


ggplot(data=data2, aes(x=stage, y=size)) +
  geom_boxplot(binaxis='y', stackdir='center') +
  theme(axis.text.x = element_text(angle=90)) 

#stages are not well defined... 

#visualize date-size relationship
#need to convert absolute dates to free time (e.g. days)

#just for grandiflorum_1_1_1####

data3<-
  read.csv('epimedium_growth_data_pivot.csv', header=TRUE) %>%
  filter(grepl('grandiflorum_1_1_3', Species_Individual_Panicle_Flower)) %>%
  na.omit() #change species name to generate separate graphs

qplot(data=data3, x=days, y=size)# +
geom_dotplot(dotsize=0.7, binwidth = 0.) +
  ylim(0,30) +
  ggtitle('Epimedium grandiflorum, flower 1')

qplot(x=x, y=y, data=test) 