  # sequence assembly: 
# for creating long consensus sequence from short fragments of same sequence

/Linux/velvet/1.2.10/velveth grandiflorum_asm 1 -fasta -short gran_seqs.fasta 

/Linux/velvet/1.2.10/velvetg grandiflorum_asm -exp_cov auto -cov_cutoff 2 -min_contig_lgth 2

-min_contig_lgth 1
-alignments yes
-min_pair_count 2

