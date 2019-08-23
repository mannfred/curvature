# #alignment seeks to perform the act
# of taking multiple divergent biological sequences of the same "type" and
# fitting them to a form that reflect some shared quality. Wright 2019


#FASTASTUFF####

data6<-
  data5 %>% 
  filter(Species_epithet == 'grandiflorum') %>%
  group_by(Species_Individual_Panicle_Flower) %>%
  pivot_wider(Species_Individual_Panicle_Flower, 
              names_from=date, 
              values_from = new_stage) %>% 
  replace(., is.na(.), "") %>%
  unite(seq, 2:14, sep="", remove=FALSE)

#write a fasta file from tibble
Xfasta <- character(nrow(data6) * 2) #empty character vector with slots for fasta header and accompanying seq
Xfasta[c(TRUE, FALSE)] <- paste0(">", data6$Species_Individual_Panicle_Flower) #paste in IDs
Xfasta[c(FALSE, TRUE)] <- data6$seq #paste in seq data

writeLines(Xfasta, "gran_seqs.fasta")  

#multiple sequence alignment####

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
names(grandiflorum_stringset) = paste(data5 %>% 
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
