#!/bin/bash

cd ~/Documents/GitHub/Bioinformatics/Midterm2/modeltest-ng

./bin/modeltest-ng \
	-i ~/Documents/GitHub/Bioinformatics/Midterm2/metazoa_alignment.5k.fasta \
	-d nt \
	-T raxml \
	-p 4 \
	-o metazoa_modeltest