############################################################
# this R script includes initial analyses of stage-size data 
#collected from Epimedium spp., April 2019, UBC Botanical Garden

library(emmeans)
library(here)
library(lme4)
library(lmerTest)
library(lubridate)
library(stringr) 
library(tidyverse)
here = here::here #mask lubridate::here




#########################################
# power analysis for epimedium experiment

pwr.t.test(d=0.30, sig.level = 0.05, power=0.7, type="two.sample", alternative = 'greater')

#Two-sample t test power calculation 

#n = 105.253, d = 0.3, sig.level = 0.05, power = 0.7, alternative = greater
#n is number of samples needed per species




############################################
# Epimedium growth data collected April 2019
# pivoting and plotting garden data 


#pivot original data frame so that "stage" info can be manually inputted 
data <-
  read.csv(
      here("data/epimedium_growth_data.csv"), 
      header=TRUE) %>%
  select(Species_Individual_Panicle_Flower, 2:25) %>% #add ", 2:25" parameter to select data upto column 25 (col25=May 2, 2019)
  as_tibble() %>%
  gather(key="date", value="size", -Species_Individual_Panicle_Flower) #pivot


levels <- 
  unique(data$Species_Individual_Panicle_Flower) %>%
  str_sort(., numeric = TRUE) #create a levels vector with the identifiers in the right order (stringr)


data_sort <-
  data %>%
  mutate(Species_Individual_Panicle_Flower = factor(Species_Individual_Panicle_Flower, levels=levels)) %>% #groups by SIPF
  arrange(Species_Individual_Panicle_Flower) %>% #arranges observations by individual
  print(n=50) #%>%
  #write.csv(
  #   here("data/epimedium_growth_data_pivot.csv"), 
  #   row.names=F) #stage info now manually entered into this file





###########################################
#visualize and test stage-size relationship

#tidying data
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
  ungroup() 
                              

#isolate numeric identifiers
data2ids <-
  data2 %>% 
  filter(Species_epithet == "koreanum" | Species_epithet == "violaceum") %>% 
  mutate(identity = str_replace(Species_Individual_Panicle_Flower,  "[a-z]+", "")) %>% #removes alphabet and leaves identifiers 
  mutate(identity = str_replace(identity, substr(identity, 1, 1), ""))  #remove leading "_" from ID strings

#add identifiers to a list
ID_matrix <- str_split_fixed(data2ids$identity, "_", n=3) #matrix of IDs in SIPC format

#add identifiers back to tibble 
data2ids <-
  data2ids %>%
  mutate(flower_ID = ID_matrix[,3]) %>% #create column for flower ID
  mutate(panicle_ID = ID_matrix[,2]) %>% #create column for panicle ID
  mutate(indiv_ID = ID_matrix[,1]) 
  




####################################################
#test for differences between stages (within species)

size_stage_model <- 
  lmerTest::lmer(size ~ stage + Species_epithet + (1|indiv_ID) , data = data2ids)

emmeans(size_stage_model, list(pairwise ~ stage + Species_epithet), adjust = "tukey")
#some stages are not distinct (within species)
#interpretation of row 1 of emmeans() output:
#the mean size of all E. koreanum at stage "E" is 2.32 mm.







############################################################
#create new stage definitions (i.e. collapse non-sig stages) 

data5<-
  data2ids %>%
  ungroup() %>%
  mutate(date = data2ids$date %>%
           paste(., "_19", sep= "") %>% #add year 2019 to date format
           mdy()) %>%
  mutate(new_stage = 
           case_when(c(Species_epithet == 'koreanum' | Species_epithet == 'violaceum') & stage == "E" ~ "C",
                     c(Species_epithet == 'koreanum' | Species_epithet == 'violaceum') & stage == "G" ~ "G",
                     stage == "P" | stage == "D" | stage == "A" ~ "A",
                     stage == "O" ~ "T" ))

#save this frequently used tibble as an .rds file
#saveRDS(data5, file="epimedium_growth_data_pivot_redefined_stages.rds") 






###########################################################
#test for size differences between "new" stages (between species)


model1 <- 
  lmerTest::lmer(size ~ Species_epithet + new_stage + (1|indiv_ID), data = data5)

qqnorm(resid(model1))
qqline(resid(model1))

tukey_results<-
  emmeans(model1, list(pairwise ~ Species_epithet + new_stage), adjust = "tukey")
 
#stages do not differ in size between species


#plot model

emmip(model1, Species_epithet~new_stage, type="response") +
  geom_point(
    aes(x = new_stage, y = size, colour = Species_epithet), 
    data = data5, 
    pch = 20, 
    alpha=0.35, 
    size=3, 
    position=position_jitterdodge(dodge.width=0.3)) +
  labs(y = "Sepal size (mm)",
       x = "Developmental stage") +
  scale_x_discrete(limits=c("C", "G", "T", "A")) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.20)) +
  stat_summary(
    fun.y = mean, 
    geom = "errorbar", 
    aes(ymax = ..y.., ymin = ..y.., group = factor(Species_epithet)),
    width = 0.5, 
    linetype = "solid", 
    position = position_dodge(),
    size=1.5,
    alpha=0.65)
  
 
#for boxplot add..
# geom_boxplot(
#   aes(x=new_stage, y=size, colour=Species_epithet), data=data5)
# theme_classic() 



#################################################
#plot data with newly defined stages using ggplot

# ggplot(
#     data=data5 %>% filter(Species_epithet=="koreanum" | Species_epithet=="violaceum"), #remove filter to include grandiflorum
#     aes(x=new_stage, y=size)) +
#     geom_boxplot(
#       aes(fill = factor(Species_epithet), #colour by species
#           factor(new_stage, levels=c("C", "G", "T", "A")))) + #reorder
#     stat_summary(
#       fun.y = mean, 
#       geom = "errorbar", 
#       aes(ymax = ..y.., ymin = ..y.., group = factor(Species_epithet)),
#       width = 0.75, 
#       linetype = "dashed", 
#       position = position_dodge()) +
#     theme_classic()  #removes gray backdrop













