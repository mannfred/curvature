#power analysis for epimedium experiment####

pwr.t.test(d=0.30, sig.level = 0.05, power=0.7, type="two.sample", alternative = 'greater')

#Two-sample t test power calculation 

#n = 105.253 
#d = 0.3
#sig.level = 0.05
#power = 0.7
#alternative = greater

#NOTE: n is number of samples needed per species




#plotting garden data####

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
                              
  

#plot sizes at varying stages
ggplot(data=data2,
         aes(
           x=stage, 
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
#E = B stage ('Bud Burst')
#G = G stage ('growth stage') 
#O = O stage ('opening stage')
#P-D-A = A stage ('anthesis')


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
                               ~ "B",
                               c(Species_epithet == 'koreanum' |
                               Species_epithet == 'violaceum') &
                               stage == "G"
                               ~ "G",
                               stage == "P" |
                               stage == "D" |
                               stage == "A"
                               ~ "A",
                               stage == "O"
                               ~ "O"
                               ))



#how does size develop through time?####

data5<-
  data2 %>%
  filter(grepl('grandiflorum_1_1_3', Species_Individual_Panicle_Flower)) %>%
  na.omit() #change species name to generate separate graphs

qplot(data=data5, x=days, y=size)# +
geom_dotplot(dotsize=0.7, binwidth = 0.) +
  ylim(0,30) +
  ggtitle('Epimedium grandiflorum, flower 1')


  
  
#dynamic time warping! http://marcocuturi.net/GA.html

dtwdata<-
  read.csv(
    here("data/test.csv"), 
    header=FALSE)

align<-dtw(dtwdata, slice(dtwdata, 3), open.end=TRUE, open.begin=TRUE, step.pattern=asymmetric, keep=TRUE)
str(align)

align$jmin

dtw_avg <- DBA(CharTraj[1:5], CharTraj[[1]], trace = TRUE)



#multiple sequence alignment
grandiflorum_stringset<-
  data5 %>% 
  filter(Species_epithet == 'grandiflorum') %>%
  group_by(Species_Individual_Panicle_Flower) %>%
  pivot_wider(Species_Individual_Panicle_Flower, 
              names_from=date, 
              values_from = new_stage) %>% 
  replace(., is.na(.), "") %>%
  unite(seq, 2:14, sep="", remove=FALSE) %>% #unites nucleotides from columns 2:14 into a single seq
  pull(seq) %>% #isolate the seqs column
  AAStringSet() #creats AAStringSet for msa

#add names
names(grandiflorum_stringset) = paste(data5 %>% 
                                        filter(Species_epithet == 'grandiflorum') %>%
                                        pull(Species_Individual_Panicle_Flower) %>%
                                        unique(.),
                                        
                                      sep=""
                                      )


library(msa) 

#for Epimedium grandiflorum
align<-msa::msaClustalW(grandiflorum_stringset, 
                        cluster="nj",
                        maxiters = 1000,
                        gapOpening = 100, #terminal gaps are not penalized
                        gapExtension = 20, 
                        type="protein")

detail(align) #inspect
str(align) #S4 'formal class MsaAAMultipleAlignment' w 6 slots




#compute consensus seq for the two (obvious) groups (1:39 and 40:46) and then 
#re-align the two consensus seqs


#convert msa object (S4) to a tibble for splitting
aligntb<-
  align %>%
  as.matrix() %>%
  as_tibble %>%
  mutate(ID = align@unmasked@ranges@NAMES) %>%
  select(ID, everything()) #move ID column to the beginning of the tibble
   
#seq group 1
align1<-
  aligntb %>%
  dplyr::slice(1:39) %>%
  unite(seq, 2:18, sep="", remove=FALSE) %>%
  pull(seq) %>% #isolate the seqs column
  AAStringSet()  #convert to XString object for msa

names(align1) = paste(aligntb %>% 
                        pull(ID) %>% 
                        .[1:39], 
                      sep=""
                      )

align1.2<-msa::msaClustalW(align1, 
                        cluster="nj",
                        maxiters = 1000,
                        gapOpening = 100, #terminal gaps are not penalized
                        gapExtension = 20, 
                        type="protein")

msaConsensusSequence(align1.2)
#[1] "----GOOOAAAAA------"

#seq group 2
align2<-
  aligntb %>%
  dplyr::slice(40:46) %>%
  unite(seq, 2:18, sep="", remove=FALSE) %>%
  pull(seq) %>% #isolate the seqs column
  AAStringSet()  #convert to XString object for msa
  
names(align2) = paste(aligntb %>% 
                        pull(ID) %>% 
                        .[40:46], 
                      sep=""
                      )

align2.2<-msa::msaClustalW(align2, 
                           cluster="nj",
                           maxiters = 1000,
                           gapOpening = 100, #terminal gaps are not penalized
                           gapExtension = 20, 
                           type="protein")

msaConsensusSequence(align2.2)
#[1] "-------GGGGGGOOAA-"



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


#write a fasta file from tibble
# fasta <- character(nrow(data6) * 2) #empty character vector with slots for fasta header and accompanying seq
# fasta[c(TRUE, FALSE)] <- paste0(">", data6$Species_Individual_Panicle_Flower) #paste in IDs
# fasta[c(FALSE, TRUE)] <- data6$seq #paste in seq data
# 
# writeLines(fasta, "gran_seqs.fasta")