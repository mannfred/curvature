# #alignment seeks to perform the act
# of taking multiple divergent biological sequences of the same "type" and
# fitting them to a form that reflect some shared quality. Wright 2019

#pairwise (overlap) and multiple sequence alignment####

library(here)
library(tidyverse)
library(msa) 

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

#[1] "--------GTTTAAAAA-"
#[1] "GGGGGGGGGTTTAAAAAA" (inferred)

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

#stage-elapsed_days boxplot 
ggplot(
  data=stages_days_sort %>% filter(stage != "-"),
  aes(x=stage, y=elapsed_days)
       ) + 
geom_boxplot(
  aes(fill = factor(stage, levels=c("G", "T", "A")))
             )+
  theme_classic() + #removes gray backdrop
  theme(legend.position="bottom")  




#pairwise overlap alignment ####
#find the longest string to use as a reference for local alignment
refseq<-
  grandiflorum_stringset[
    which.max(
      nchar(grandiflorum_stringset)
    )
    ] %>%
  as.character()
#     width seq                  names               
# [1] 12    GGGGGGGTTAAA        grandiflorum_1_2_4

pair<-pairwiseAlignment(grandiflorum_stringset, 
                        refseq, 
                        type="overlap",
                        substitutionMatrix=matchmatrix,
)

#compute consensus
conMatrix<-consensusMatrix(pair, as.prob=TRUE)[-1,] #remove row 1 to remove gaps
t(apply(conMatrix, 1, diff)) #find differences in probability 

#(G)GGGTTTTAAAA assuming state at 0 is G

conseq<-"GGGGTTTTAAAA" 

pair<-pairwiseAlignment(grandiflorum_stringset, 
                        refseq, 
                        type="overlap",
                        substitutionMatrix=matchmatrix,
)

