#!/bin/bash

#Need to run this to run this script
#chmod +x run_raxml.sh

#make directory to store raxml results
mkdir -p ~/Documents/GitHub/Bioinformatics/Midterm2/raxml_run
cd ~/Documents/GitHub/Bioinformatics/Midterm2/raxml_run

#run raxml and set arguments
raxml-ng \
	--all \
	--msa ~/Documents/GitHub/Bioinformatics/Midterm2/metazoa_alignment.5k.fasta \
	--model GTR+I+G4 \
	--prefix metazoa \
	--threads 4 \
	--seed 42 \
	--bs-trees 100

# arguments
#all = all-in-one (ML search + bootstrapping)
#msa = alignment file
#model = model specifications (decided by modeltest)
#prefix = prefix for output file name
#threads = number of parallel threads to use
#seed = pseudo-random number generator (42 is the answer to everything)
#bs-trees = number of bootstrap replicates (100 is standard minimum)

#raxml was chosen due to support from ModelTest-NG.
#It ran maximum likelihood to generate a single best tree with bootstrapping.
