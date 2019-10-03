library(here)
library(tidyverse)
library(msa) 

koreanum_stringset<-
  data5 %>% 
  filter(Species_epithet == 'koreanum') %>%
  group_by(Species_Individual_Panicle_Flower) %>%
  pivot_wider(Species_Individual_Panicle_Flower, 
              names_from=date, 
              values_from = new_stage) %>% 
  replace(., is.na(.), "") %>%
  unite(seq, 2:14, sep="", remove=FALSE) %>% #unites nucleotides from columns 2:14 into a single seq
  pull(seq) %>% #isolate the seqs column
  AAStringSet() #creats AAStringSet for msa


#add names to stringset
names(koreanum_stringset) = paste(
  data5 %>% 
    filter(Species_epithet == 'koreanum') %>%
    pull(Species_Individual_Panicle_Flower) %>%
    unique(.),
  sep=""
                                  )

matchmatrix<-
  read.table(here("data/match_matrix.txt")) %>% 
  as.matrix 

colnames(matchmatrix)[24]<-"*" #gets turned into ".X" during read.table for some reason..

#multiple sequence alignment
align<-msa::msaClustalW(koreanum_stringset, 
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

#[1] "------------GG--T-AA----"
#[1] "CCCCCCCCCCCCGGGTTTAAAAAA" (inferred from consensusMatrix(align) by disregarding gaps)

#merge alignment with days 
stages_days<-
  align %>%
  as.matrix() %>%
  as_tibble() %>%
  rename_at(colnames(.[,1:24]), ~ as.character(c(1:24))) %>% #replace V1:V24 with 1:24 to be used as "number of days"
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


