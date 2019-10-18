#with stages now defined, censored and fragmented development data
#needs to be aligned to determine the "average" developmental sequence of events


# #alignment seeks to perform the act
# of taking multiple divergent biological sequences of the same "type" and
# fitting them to a form that reflect some shared quality. Wright 2019


#######################################################
#pairwise (overlap) and multiple sequence alignment####

library(here)
library(tidyverse)
library(msa) 

grandiflorum_stringset<-
  data5 %>% 
  filter(Species_epithet == 'grandiflorum') %>%
  group_by(Species_Individual_Panicle_Flower) %>%
  pivot_wider(Species_Individual_Panicle_Flower, #pivots developmental sequence data horizontally (but each step is in a separate cell)
              names_from=date, 
              values_from = new_stage) %>% 
  replace(., is.na(.), "") %>%
  unite(seq, 2:14, sep="", remove=FALSE) %>% #unites stages from columns 2:14 into a single seq
  pull(seq) %>% #isolate the seqs column
  AAStringSet() #creats AAStringSet for msa (from character vectors)


#add names to stringset
names(grandiflorum_stringset) = paste(
                                  data5 %>% 
                                  filter(Species_epithet == 'grandiflorum') %>%
                                  pull(Species_Individual_Panicle_Flower) %>%
                                  unique(.),
                                  sep=""
                                      )


#identity substiution matrix from NCBI  ftp://ftp.ncbi.nih.gov/blast/matrices/
matchmatrix<-
  read.table(here("data/match_matrix.txt")) %>% 
  as.matrix 

colnames(matchmatrix)[24]<-"*" #gets turned into ".X" during read.table for some reason..



#multiple sequence alignment
align<-msa::msaClustalW(grandiflorum_stringset, 
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

consensusMatrix(align)

#[1] "--------GTTTAAAAA-"
#[1] "GGGGGGGGGTTTAAAAAA" (inferred from consensusMatrix(align) by disregarding gaps)


##################################
#merge alignment with elapsed time

stages_days<-
  align %>%
  as.matrix() %>%
  as_tibble() %>%
  rename_at(colnames(.[,1:18]), ~ as.character(c(1:18))) %>% #replace V1:V18 with 1:18 to be used as "number of days"
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


############################################
#stage-elapsed_days boxplot (5 value spread)

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

saveRDS(stages_days_sort, file="stages_days_sort_grandiflorum.rds")



#######################################################
#calculate mean + variance for stages 
#on average, how many days have elapsed for each stage?

stages_days_sort <- readRDS("stages_days_sort_grandiflorum.rds")


#mean elapsed days per stage

library(gmodels) #for CI

stages_days_sort %>%
  group_by(stage) %>%
  na_if("-") %>%
  summarise(mean = ci(elapsed_days)[1],
            loCI = ci(elapsed_days)[2],
            hiCI = ci(elapsed_days)[3],
            stdv = ci(elapsed_days)[4])

#G: 7.01 +/- 0.45 days
#T: 10.9 +/- 0.20 days
#A: 14.7 +/- 0.20 days


lm_test<- 
  lm(elapsed_days ~ stage, data = stages_days_sort %>% na_if("-"))
     
TukeyHSD(aov(lm_test))
# elapsed days per stage are sig different (p=0)

# elapsed_days-stage scatterplot
# ggplot(
#   data=stages_days_sort %>% filter(stage != "-"), #remove gaps "-"
#   aes(x=elapsed_days, y=stage)
# ) +
# geom_point(
#   aes(x=elapsed_days, y=factor(stage, levels=c("C", "G", "T", "A")))
#            )
#   

