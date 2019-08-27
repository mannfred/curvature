#power analysis for epimedium experiment####

pwr.t.test(d=0.30, sig.level = 0.05, power=0.7, type="two.sample", alternative = 'greater')

#Two-sample t test power calculation 

#n = 105.253, d = 0.3, sig.level = 0.05, power = 0.7, alternative = greater
#n is number of samples needed per species




#rearranging and plotting garden data####

library(stringr) 
library(tidyverse)
library(lubridate)
library(here)
here = here::here #mask lubridate::here

#pivot original data frame so that "stage" info can be manually inputted 

data<-
  read.csv(
      here("data/epimedium_growth_data.csv"), 
      header=TRUE) %>%
  select(Species_Individual_Panicle_Flower, 2:25) %>% #add ", 2:25" parameter to select data upto column 25 (col25=May 2, 2019)
  as.tibble() %>%
  gather(key="date", value="size", -Species_Individual_Panicle_Flower) #pivot

levels<- 
  unique(data$Species_Individual_Panicle_Flower) %>%
  str_sort(., numeric = TRUE) #create a levels vector with the identifiers in the right order (stringr)

data_sort<-
  data %>%
  mutate(Species_Individual_Panicle_Flower = factor(Species_Individual_Panicle_Flower, levels=levels)) %>%
  arrange(Species_Individual_Panicle_Flower) %>% #arranges observations by individual
  print(n=50) %>%
  write.csv(
      here("data/epimedium_growth_data_pivot.csv"), 
      row.names=F) #stage info now manually entered into this file


#visualize and test stage-size relationship####

#data wrangling
data2<-
  read.csv(
      here("data/epimedium_growth_data_pivot.csv"), 
      header=TRUE) %>%
  na.omit() %>%
  mutate(stage = factor(stage, levels=c("E", "G", "O", "P", "D", "A"))) %>% #reorder $stage
  mutate(
    Species_epithet = str_extract(
                         Species_Individual_Panicle_Flower, 
                         "[a-z]+") #extacts words from strings
        ) %>%
  group_by(Species_Individual_Panicle_Flower) %>% #creates groups by individual 
  mutate(days = row_number()) %>% #for every level of SIPF, assign "day" numbers starting from 1
  ungroup() %>% #because I don't actually want them grouped by individual, I just did this for generating the "days" column
  group_by(Species_epithet) #grouping by grandiflorum, koreanum, violaceum
                              
  

#test for differences between stages (within species)

#lm function for map()
lm_size_stage<- 
  function(df) {
    lm(size ~ stage, data = df)
                }

#tukey function for map()
tukey_func<-
  function(lm_fit) {
    TukeyHSD(aov(lm_fit))
                    }
data3<-
  data2 %>%
  group_by(Species_epithet) %>%
  nest() %>% #nests by Species_epithet
  mutate(lm_fit = 
           map(data, lm_size_stage) #fits lms 
         ) %>%
  mutate(tukey =
           map(lm_fit, tukey_func) #compares means within species
         ) 




#test for size differences between stages (between species)

#lm function for map()
lm_size_species<- 
  function(df) {
    lm(size ~ Species_epithet, data = df)
                }

data4<-
  data2 %>%
  group_by(stage) %>% #grouping by stage to make between-species comparisons
  nest() %>%
  mutate(lm_fit = 
           map(data, lm_size_species)
         ) %>%
  mutate(tukey =
           map(lm_fit, tukey_func)
         )



#create new stage definitions (i.e. collapse non-sig stages) 
#E = C stage ('child stage')
#G = G stage ('growth stage') 
#O = T stage ('nectar producTion stage')
#P-D-A = A stage ('anthesis')
#use https://www.hiv.lanl.gov/content/sequence/HelpDocs/IUPAC.html for checking AA codes

data5<-
  data2 %>%
  ungroup() %>%
  mutate(date = data2$date %>%
           paste(., "_19", sep= "") %>% #add year 2019 to date format
           mdy() 
         ) %>%
  mutate(new_stage = case_when(c(stage == "E" | 
                               stage == "G") & 
                               Species_epithet == 'grandiflorum' 
                               ~ "G", 
                               c(Species_epithet == 'koreanum' |
                               Species_epithet == 'violaceum') &
                               stage == "E" 
                               ~ "C",
                               c(Species_epithet == 'koreanum' |
                               Species_epithet == 'violaceum') &
                               stage == "G"
                               ~ "G",
                               stage == "P" |
                               stage == "D" |
                               stage == "A"
                               ~ "A",
                               stage == "O"
                               ~ "T"
                               ))


#plot sizes at varying stages
ggplot(data=data5,
       aes(
         x=new_stage, 
         y=size 
       )
) +
  geom_boxplot(
    aes(
      fill = factor(Species_epithet)
    )
  )+
  theme(axis.text.x = element_text(angle=90)) +
  theme_classic() + #removes gray backdrop
  theme(legend.position="bottom")  












ggplot(data=data6,
       aes(
         x=date, 
         y=size, 
         group = Species_Individual_Panicle_Flower,
         colour = factor(Species_epithet)
       )
) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90)) +
  theme_classic() + #removes gray backdrop
  theme(legend.position="bottom") 
#ordinal date



#how does size develop through time?####
# 
# data5<-
#   data2 %>%
#   filter(grepl('grandiflorum_1_1_3', Species_Individual_Panicle_Flower)) %>%
#   na.omit() #change species name to generate separate graphs
# 
# qplot(data=data5, x=days, y=size)# +
# geom_dotplot(dotsize=0.7, binwidth = 0.) +
#   ylim(0,30) +
#   ggtitle('Epimedium grandiflorum, flower 1')


