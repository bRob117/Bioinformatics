#Bryce Robinson
#BSC4434C Bioinformatics
#Lab 10 Protein
#Using ChatGPT when stuck with no generated comments

#set working directory
setwd("~/FGCU/BSC4434C/Lab10")

#install new protein packages
#install.packages("UniprotR")
#install.packages("protti")
#BiocManager::install("GenomicAlignments")
#install.packages("r3dmol")

#load packages
library(UniprotR)
library(protti)
library(r3dmol)
library(GenomicAlignments)
library(Biostrings)

#load amino acid sequence for lycopene beta cyclase gene as string set
AAseq <- readAAStringSet("lycopeneAA.fasta")

#Solanum lycopersicum lycopene beta cyclase gene UniProt five accession matches (actual accession for gene is Q43503)
#A0ABM1GMN2
#A0AAF0QJQ2
#A0A9J5ZDD7
#A0ABD2UT19
#A0AAV9L8J8

#manually wrote .csv from UniProt results
#read text file while ignoring first row information (kept getting error about empty columns)
accessions <- read.csv("lycopene_accessions.csv", comment.char = "#", header = FALSE)

#convert all variables to a single vector (vector needed for GetProteinGOInfo function)
accessions_character <- as.character(unlist(accessions))

#convert all values into one string (assignment said to convert into string, idk why)
accessions_string <- paste(accessions_character, collapse = ",")

#get Gene Ontology terms and print
go_terms <- GetProteinGOInfo(accessions_character)
print(go_terms)
#pasted results:
#Gene.Ontology..molecular.function.
#A0ABM1GMN2                                                                                                                                                          <NA>
#A0AAF0QJQ2 neoxanthin synthase activity [GO:0034020]; oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen [GO:0016705]
#A0A9J5ZDD7 neoxanthin synthase activity [GO:0034020]; oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen [GO:0016705]
#A0ABD2UT19                                                                                                                     neoxanthin synthase activity [GO:0034020]
#A0AAV9L8J8 neoxanthin synthase activity [GO:0034020]; oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen [GO:0016705]
#Gene.Ontology..cellular.component.
#A0ABM1GMN2                               <NA>
#A0AAF0QJQ2           chloroplast [GO:0009507]
#A0A9J5ZDD7           chloroplast [GO:0009507]
#A0ABD2UT19           chloroplast [GO:0009507]
#A0AAV9L8J8           chloroplast [GO:0009507]

#plot results
PlotGoInfo(go_terms)

#Handy visualization for publications code taken from UniprotR GitHub (idk what is supposed to be in my plot exactly so not sure if it is right)
PlotGOAll(GOObj = go_terms, Top = 10, directorypath = getwd(), width = 8, height = 5)

#use provided functions to find info on diseases or pathologies
go_pathology <- GetPathology_Biotech(accessions_character)
go_diseases <- Get.diseases(accessions_character)
print(go_pathology)
print(go_diseases)
#results are NA across the board. most likely since is a plant gene and isn't associated with human diseases.

#access structural info using fetch function in protti package 
protein_structure <- fetch_uniprot(accessions_character)

#extract pdb out of fetch
pdb_ids <- protein_structure$xref_pdb
#NA for all five pdbs

#pull structural info from Protein DataBase for "1ZMR" and "2HWG" since no data were available
ZMR <- fetch_pdb("1ZMR")
HWG <- fetch_pdb("2HWG")

#get info on 3D structure and visualize using alphafold website
alphafold <- fetch_alphafold_prediction("Q43503")

#this assignment would have been better had I chose a better gene but I really don't want to redo the whole thing. It is cool seeing the predicted protein folded.
