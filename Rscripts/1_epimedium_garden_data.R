# this R script includes initial analyses of stage-size data 
#collected from Epimedium spp., April 2019, UBC Botanical Garden



############################################
#power analysis for epimedium experiment####

pwr.t.test(d=0.30, sig.level = 0.05, power=0.7, type="two.sample", alternative = 'greater')

#Two-sample t test power calculation 

#n = 105.253, d = 0.3, sig.level = 0.05, power = 0.7, alternative = greater
#n is number of samples needed per species




############################################
# Epimedium growth data collected April 2019
# rearranging and plotting garden data ##### 

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
  as_tibble() %>%
  gather(key="date", value="size", -Species_Individual_Panicle_Flower) #pivot

levels<- 
  unique(data$Species_Individual_Panicle_Flower) %>%
  str_sort(., numeric = TRUE) #create a levels vector with the identifiers in the right order (stringr)

data_sort<-
  data %>%
  mutate(Species_Individual_Panicle_Flower = factor(Species_Individual_Panicle_Flower, levels=levels)) %>% #groups by SIPF
  arrange(Species_Individual_Panicle_Flower) %>% #arranges observations by individual
  print(n=50) #%>%
  #write.csv(
  #   here("data/epimedium_growth_data_pivot.csv"), 
  #   row.names=F) #stage info now manually entered into this file






###############################################
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

#creating functions to pass onto purrr::map()
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

#creating tibble with comparison of mean sizes between stages (within species)
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
#creating tibble with comparison of mean sizes between stages (between species)
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

data5<-
  data2 %>%
  ungroup() %>%
  mutate(date = data2$date %>%
           paste(., "_19", sep= "") %>% #add year 2019 to date format
           mdy() 
         ) %>%
  mutate(new_stage = case_when(c(stage == "E" |                  #for E. grandiflorum, collapse E and G in to G
                               stage == "G") & 
                               Species_epithet == 'grandiflorum' 
                               ~ "G", 
                               c(Species_epithet == 'koreanum' | #for E. kore and E. viol, convert E to C
                               Species_epithet == 'violaceum') &
                               stage == "E" 
                               ~ "C",
                               c(Species_epithet == 'koreanum' | #keep G stage for E. kore and E. viol
                               Species_epithet == 'violaceum') &
                               stage == "G"
                               ~ "G",
                               stage == "P" |                    #collapse PDA into A 
                               stage == "D" |
                               stage == "A"
                               ~ "A",
                               stage == "O"                      #convert O to T
                               ~ "T"
                               )) 

#save this frequently used tibble as an .rds file
saveRDS(data5, file="epimedium_growth_data_pivot_redefined_stages.rds") 



####################################
#plot data with newly defined stages
#plot sizes at varying stages

ggplot(
  data=data5 %>% filter(Species_epithet=="koreanum" | Species_epithet=="violaceum"), #remove filter to include grandiflorum
  aes(x=new_stage, y=size)
       ) +
geom_boxplot(
  aes(fill = factor(Species_epithet), #colour by species
      factor(new_stage, levels=c("C", "G", "T", "A")) #reorder
      )
             ) +
theme(axis.text.x = element_text(angle=90)) +
theme_classic() + #removes gray backdrop
theme(legend.position="bottom")  




##########################################################
#rerun comparison of mean sizes with new stage definitions


#lm function for map()
lm_size_new_stage<- 
  function(df) {
    lm(size ~ new_stage, data = df)
  }


#tukey function for map()
tukey_func<-
  function(lm_fit) {
    TukeyHSD(aov(lm_fit))
  }

#tibble with model fits and tukey's HSD results
data6<-
  data5 %>%
  group_by(Species_epithet) %>%
  nest() %>% #nests by Species_epithet
  mutate(lm_fit = 
           map(data, lm_size_new_stage) #fits lms 
  ) %>%
  mutate(tukey =
           map(lm_fit, tukey_func) #compares means within species
  ) 


# sepal size per stage are sig different (p=0 for all within species comparisons). 





############################################
# calculate mean +/- CI sizes for all 3 spp.
# what is the size of the flower at each stage?

library(gmodels) #for CI

data5_size_summary<-
  readRDS(here("/data/RDS_files/epimedium_growth_data_pivot_redefined_stages.rds")) %>%
  group_by(Species_epithet, new_stage) %>% #group by species, then group by new_stage
  summarise(mean = ci(size)[1],            #ci() function produces 4 results, here we assign them labels
            loCI = ci(size)[2],
            hiCI = ci(size)[3],
            stdv = ci(size)[4])

write.csv(data5_size_summary, 
          row.names=FALSE,
          file=here("/data/Table_4_size_stage_summary.csv"))




##########################################
# testing for independent sampling units #

library(tidyverse)
library(stringr)
library(here)

#isolate numeric identifiers
data5<-
  readRDS(here("/data/RDS_files/epimedium_growth_data_pivot_redefined_stages.rds")) %>%
  mutate(identity = str_replace(Species_Individual_Panicle_Flower,  "[a-z]+", "")) %>% #removes alphabet and leaves identifiers 
  mutate(identity = str_replace(identity, substr(identity, 1, 1), ""))  #remove leading "_" from ID strings
  
#add identifiers to a list
ID_matrix<-str_split_fixed(data5$identity, "_", n=3) #matrix of IDs in SIPC format


data5_ids<-
  data5 %>%
  mutate(flower_ID = ID_matrix[,3]) %>% #create column for flower ID
  mutate(panicle_ID = ID_matrix[,2]) %>% #create column for panicle ID
  mutate(indiv_ID = ID_matrix[,1]) 




size_stage_model <- lmer(size ~ days + (1|indiv_ID), data=data5_ids %>% filter(Species_epithet == 'koreanum'))



