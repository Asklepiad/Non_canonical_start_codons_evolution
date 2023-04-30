#!/bin/bash

input=$1
threads=$2
curr_dir="../${input}/data/evolution_models"
if [[ ! -e ${curr_dir} ]]
then
        mkdir ${curr_dir};
fi

directory="multialignments/"


for f in $(ls ../${input}/data/${directory}*.fasta); do
	if [[ $length -lt 1000 ]]
        then
                muscle -align $f -threads $threads -output ${f%.fasta}.afa;
        else
                muscle -super5 $f -threads $threads -output ${f%.fasta}.afa;
        fi
done

