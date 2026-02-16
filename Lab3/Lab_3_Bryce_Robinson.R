#set working directory
setwd("~/FGCU/BSC4434C/Lab3")

#load packages
library(Biostrings)
library(pwalign)

#import DNA sequences
Seq20 <- readDNAStringSet("Sequence20")
SeqG_hirsutus <- readDNAStringSet("G_hirsutus_mitochondria")

#align sequences
alignment_results <- pairwiseAlignment(Seq20, SeqG_hirsutus, type = "local")

#view alignment
print(alignment_results)
