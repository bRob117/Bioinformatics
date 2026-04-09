#!/bin/bash

mkdir -p ~/Documents/GitHub/Bioinformatics/Midterm2/raxml_run
cd ~/Documents/GitHub/Bioinformatics/Midterm2/raxml_run

raxml-ng \
	--all \
	--msa ~/Documents/GitHub/Bioinformatics/Midterm2/metazoa_alignment.5k.fasta \
	--model GTR+I+G4 \
	--prefix metazoa \
	--threads 4 \
	--seed 42 \
	--bs-trees 100