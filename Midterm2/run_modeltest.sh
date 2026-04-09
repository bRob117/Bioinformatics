#!/bin/bash

#Need to run this to run this script
#chmod +x run_modeltest.sh

#change directory
cd ~/Documents/GitHub/Bioinformatics/Midterm2/modeltest-ng

#run modeltest and set arguments
./bin/modeltest-ng \
	-i ~/Documents/GitHub/Bioinformatics/Midterm2/metazoa_alignment.5k.fasta \
	-d nt \
	-T raxml \
	-p 4 \
	-o metazoa_modeltest

#arguments
#i = input alignment file
#d = datatype -> nt = nucleotide
#T = template -> raxml = RAxML
#p = threads -> 4
#o = output file