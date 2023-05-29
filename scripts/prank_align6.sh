#!/bin/bash


# Reads arguments
input=$1


# Creates folder for storing data about best evolution model for phylogenetic trees
curr_dir="../${input}/data/evolution_models"
if [[ ! -e ${curr_dir} ]]
then
        mkdir ${curr_dir};
fi


# Multialigns ortologous rows with inconsistent start codon
directory="multialignments/"
for f in $(ls ../${input}/data/${directory}*.fasta); do
        prank -d=$f -o=${f%.fasta}_prank;
done
