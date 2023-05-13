#!/bin/bash

input=$1
curr_dir="../${input}/data/evolution_models"
if [[ ! -e ${curr_dir} ]]
then
        mkdir ${curr_dir};
fi

directory="multialignments/"
for f in $(ls ../${input}/data/${directory}*.fasta); do
        prank -d=$f -o=${f%.fasta}_prank;
done
