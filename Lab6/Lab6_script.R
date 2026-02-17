#set working directory
setwd("~/FGCU/BSC4434C/GitHub/Bioinformatics/Lab6")

#load multiple sequence alignment package
library(msa)
library(seqinr)

#import sequence files
e_coli_seqs <- readDNAStringSet("seqdump.fasta")

#visualize sequences
print(e_coli_seqs)

#align sequences with MUSCLE
aligned_seqs <- msa(e_coli_seqs, method="Muscle")

#check for gaps
print(aligned_seqs, show="complete")

#measure consensus length
e_coli_consensus <- consensusString(aligned_seqs)
nchar(e_coli_consensus)

#calculate GC content
consensus_vector <- s2c(e_coli_consensus)
GC(consensus_vector)

#compute distance matrix
dist_matrix <- consensusMatrix()
