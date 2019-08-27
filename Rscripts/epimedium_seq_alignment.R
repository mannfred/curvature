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

#for visualization of the alignment process
align<-msa::msaClustalW(grandiflorum_stringset, 
                        cluster="nj", #neighbour joining 
                        maxiters = 1000,
                        gapOpening = 100, #terminal gaps are not penalized
                        gapExtension = 20, 
                        substitutionMatrix = matchmatrix,
                        type="protein")

detail(align) #visualize


?msaConsensusSequence
msaConsensusSequence(align, 
                     type="upperlower", 
                     thresh=c(10, 0.0000000001), 
                     ignoreGaps=FALSE)

#[1] "--------GTTTAAAAA-"
#[1] "GGGGGGGGGTTTAAAAAA" (inferred)
