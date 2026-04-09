setwd("~/Documents/GitHub/Bioinformatics/Midterm2")

library(ape)
library(ggtree)
library(treeio)

tree <- read.tree("raxml_run/metazoa.raxml.support")

rooted_tree <- ape::root(tree, outgroup = c("Plakina_jani", "Grantia_compressa"), resolve.root = TRUE)

plot(rooted_tree, show.node.label = TRUE)

p <- ggtree(rooted_tree) +
  geom_tiplab(fontface = "italic", size = 2.8, align = TRUE, linetype = "dotted") +
  geom_nodelab(size = 2.5, hjust = -0.2)

ggsave("metazoa_rooted_tree.pdf", p, width = 10, height = 11)

library(Biostrings)

aln <- readDNAStringSet("metazoa_alignment.gene.fasta")

human <- aln[grep("Homo_sapiens", names(aln))]

human_nogap <- DNAStringSet(gsub("-", "", as.character(human)))

human_protein <- translate(human_nogap)

writeXStringSet(human_protein, "human_protein.fasta")

library(org.Hs.eg.db)
library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)

accession <- "P54098"

go_ann <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = accession,
                                columns = c("SYMBOL","GO","ONTOLOGY"),
                                keytype = "UNIPROT")

go_ann$TERM <- AnnotationDbi::Term(GOTERM[go_ann$GO])

go_unique <- dplyr::distinct(go_ann, GO, ONTOLOGY, TERM)

go_unique |>
  dplyr::group_by(ONTOLOGY) |>
  dplyr::slice_head(n = 3) |>
  print(n = Inf)

counts <- dplyr::count(go_unique, ONTOLOGY)

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

ggsave("GO_subontology_plot.png", hp, width = 6, height = 4)