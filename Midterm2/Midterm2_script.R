#Bryce Robinson
#BSC4434C Bioinformatics
#Midterm 2 Take Home Practical for Midterm 2
#Uncommented code support from Claude 1.1348.0

#---------------------------------------------------------------------------------------------------------
#1 Estimate a phylogenetic tree from aligned cDNA using ModelTest-NG and its supported model: RAxML GTR+I+G

#2 Root the tree on the branch leading to Plakina jani and Grantia compressa.

#set working directory
setwd("~/Documents/GitHub/Bioinformatics/Midterm2")

#load packages
library(ape)
library(ggtree)
library(treeio)

#load tree data into variable "tree"
tree <- read.tree("raxml_run/metazoa.raxml.support")

#set outgroups as Plakina jani and Grantia compressa using "tree" phylo object.
#c() combines values into a vector
#resolve.root specifies if new root is a bifurcating node
rooted_tree <- ape::root(tree, outgroup = c("Plakina_jani", "Grantia_compressa"), resolve.root = TRUE)

#plot tree with node labels
plot(rooted_tree, show.node.label = TRUE)

#3 Find a way to measure how “good” your tree is. Is this a good tree or a bad one? How can you tell? 
#Are there some parts of the tree that are better than others?
#Answer in Word document.

#4 There are six groups of taxa in your tree (Lophotrochozoa, Ecdysozoa, Xenambulacraria, Chordata, Cnidaria, and Porifera). 
#Does your tree support the monophyly of these clades? 
#What evidence supports or does not support your conclusion?
#Answer in Word document.

#5 Plot your rooted tree in the tree visualizer of your choice, and export a pdf or similar file (e.g., png, jpeg). 
#Include this file in your submission.

#set ggtree figure as variable "p"
#geom functions adjust sizes, positions, and fonts of tip and node labels. align makes tips equal.
#hexpand adds space on the right for tip labels
p <- ggtree(rooted_tree) +
  geom_tiplab(fontface = "italic", size = 2.8, align = TRUE, linetype = "dotted") +
  geom_nodelab(size = 2.5, hjust = -0.2)

#save figure as pdf
ggsave("metazoa_rooted_tree.pdf", p, width = 10, height = 11)

#---------------------------------------------------------------------------------------------------------
#I had to restart R to unload packages for the next half to work properly.

#6 Next, you will investigate protein function. You have been provided the data from one of the mRNAs. 
#Your goal is to investigate the function of the protein. Use the metazoa_alignment.gene.fasta file for questions 6-10.

#7.	Look at the samples in your alignment. What biological/genetic process might have resulted in some samples having so many dashes in the alignment? 
# What evidence led you to that conclusion?
#Answer in Word document.

#8.	Load your file into R. Find the gene from Homo sapiens in your alignment. Translate that sequence to protein. Write it to a .fasta file.

#load package
library(Biostrings)

#load fasta file into string set
aln <- readDNAStringSet("metazoa_alignment.gene.fasta")

#extract human sequence from others using grep
human <- aln[grep("Homo_sapiens", names(aln))]

#remove gaps from sequence using gsub
human_nogap <- DNAStringSet(gsub("-", "", as.character(human)))

#translate to amino acid sequence
human_protein <- translate(human_nogap)

#write amino acid sequence into fasta
writeXStringSet(human_protein, "human_protein.fasta")

#9 Use a database to figure out what your protein is. Click on the record for the best match. What is the name of the gene?
#Answer in Word document.

#10 Locate the closest match to your sample in the UniProt database. Find the UniProt accession number. Using an R script, find the GO sub-ontologies for this gene. 
#List at least one term from each of the three sub-ontologies. Plot the GO info and provide the plot in your submission.

#load packages
library(org.Hs.eg.db)
library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)

#set accession from UniProt for Homo sapiens as variable
accession <- "P54098"

#database query of organism.Homosapiens.entrezgene.database
go_ann <- AnnotationDbi::select(org.Hs.eg.db,
              keys = accession,
              columns = c("SYMBOL","GO","ONTOLOGY"),
              keytype = "UNIPROT")

#replace identifiers with names
go_ann$TERM <- AnnotationDbi::Term(GOTERM[go_ann$GO])

#remove duplicate annotations
go_unique <- dplyr::distinct(go_ann, GO, ONTOLOGY, TERM)

#print terms from sub-ontology
go_unique |>
  dplyr::group_by(ONTOLOGY) |>
  dplyr::slice_head(n = 3) |>
  print(n = Inf)

#count the rows in the table
counts <- dplyr::count(go_unique, ONTOLOGY)

#create plot with ONTOLOGY on the x axis and n (number of rows) on the y. geom, scale, labs, and theme are for "prettying up" the plot
hp <- ggplot(counts, aes(ONTOLOGY, n, fill = ONTOLOGY)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.3) +
  scale_x_discrete(labels = c(BP = "Biological\nProcess",
                              MF = "Molecular\nFunction",
                              CC = "Cellular\nComponent")) +
  labs(title = paste("GO sub-ontology annotations for", accession),
       x = NULL, y = "Number of GO terms") +
  theme_minimal() +
  theme(legend.position = "none")

#save the plot as a png
ggsave("GO_subontology_plot.png", hp, width = 6, height = 4)
