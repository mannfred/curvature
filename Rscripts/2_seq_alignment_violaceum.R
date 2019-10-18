library(here)
library(tidyverse)
library(msa) 


#see "2_seq_alignment_grandiflorum.R" for complete code annotations and explanations


violaceum_stringset<-
  data5 %>% 
  filter(Species_epithet == 'violaceum') %>%
  group_by(Species_Individual_Panicle_Flower) %>%
  pivot_wider(Species_Individual_Panicle_Flower, 
              names_from=date, 
              values_from = new_stage) %>% 
  replace(., is.na(.), "") %>%
  unite(seq, 2:14, sep="", remove=FALSE) %>% #unites nucleotides from columns 2:14 into a single seq
  pull(seq) %>% #isolate the seqs column
  AAStringSet() #creats AAStringSet for msa


#add names to stringset
names(violaceum_stringset) = paste(
  data5 %>% 
    filter(Species_epithet == 'violaceum') %>%
    pull(Species_Individual_Panicle_Flower) %>%
    unique(.),
  sep=""
)

matchmatrix<-
  read.table(here("data/match_matrix.txt")) %>% 
  as.matrix 

colnames(matchmatrix)[24]<-"*" #gets turned into ".X" during read.table for some reason..

#multiple sequence alignment
align<-msa::msaClustalW(violaceum_stringset, 
                        cluster="nj", #neighbour joining 
                        maxiters = 1000,
                        gapOpening = 100, #terminal gaps are not penalized
                        gapExtension = 20, 
                        substitutionMatrix = matchmatrix,
                        type="protein")

detail(align) #visualize


msaConsensusSequence(align, 
                     type="upperlower", 
                     thresh=c(10, 0.01), 
                     ignoreGaps=FALSE)

#[1] "--------CCCGGTTAA---"
#[1] "CCCCCCCCCCCGGTTAAAAA" (inferred from consensusMatrix(align) by disregarding gaps)

#merge alignment with days 
stages_days<-
  align %>%
  as.matrix() %>%
  as_tibble() %>%
  rename_at(colnames(.[,1:20]), ~ as.character(c(1:20))) %>% #replace V1:V20 with 1:20 to be used as "number of days"
  mutate(ID = align@unmasked@ranges@NAMES) %>% #create ID column
  select(ID, everything()) %>% #moves ID column to the front
  gather(key="elapsed_days", value="stage", -ID)  #pivot

levels<- 
  unique(stages_days$ID) #define individual IDs

stages_days_sort<-
  stages_days %>%
  mutate(ID = factor(ID, levels=levels)) %>% #groups by individual IDs
  arrange(ID) %>%
  mutate(elapsed_days = elapsed_days %>% as.numeric()) %>%
  print() 

saveRDS(stages_days_sort, file="stages_days_sort_violaceum.rds")

#stage-elapsed_days boxplot 
ggplot(
  data=stages_days_sort %>% filter(stage != "-"), #remove gaps "-"
  aes(x=stage, y=elapsed_days)
) + 
  geom_boxplot(
    aes(fill = factor(stage), #colour by stage
        factor(stage, levels=c("C", "G", "T", "A")) #reorder
    )
  )+
  theme_classic() + #removes gray backdrop
  theme(legend.position="bottom")   

#calculate mean + variance for stages ####
#on average, how many days have elapsed for each stage?

stages_days_sort <- readRDS("stages_days_sort_violaceum.rds")


#mean elapsed days per stage

library(gmodels) #for CI

stages_days_sort %>%
  group_by(stage) %>%
  na_if("-") %>%
  summarise(mean = ci(elapsed_days)[1],
            loCI = ci(elapsed_days)[2],
            hiCI = ci(elapsed_days)[3],
            stdv = ci(elapsed_days)[4])

#C: 8.27 +/- 0.40 days
#G: 12.0 +/- 0.20 days
#T: 14.1 +/- 0.20 days
#A: 16.7 +/- 0.30 days


lm_test<- 
  lm(elapsed_days ~ stage, data = stages_days_sort %>% na_if("-"))

TukeyHSD(aov(lm_test))
# elapsed days per stage are sig different (p=0)
