#!/bin/bash


# Reads arguments
mask=$1
radiobutton=$2
threads=$3


# Creates directory for storing phylogenetic trees
curr_dir1="../${mask}/data/raxmlng_trees"
if [[ ! -e ${curr_dir1} ]]
then
        mkdir ${curr_dir1};
fi

# Creates directory for storing config files for iTOL visualisation 
curr_dir2="../${mask}/data/configs"
if [[ ! -e ${curr_dir2} ]]
then
        mkdir ${curr_dir2};
fi


# Creates directories for the next stages of pipeline, when we wiil recompute data and will have more confindent start codons.
curr_dir3="../${mask}/data/new_multialignments"
if [[ ! -e ${curr_dir3} ]]
then
        mkdir ${curr_dir3};
fi

curr_dir4="../${mask}/data/new_orto_rows"
if [[ ! -e ${curr_dir4} ]]
then
        mkdir ${curr_dir4};
fi

curr_dir5="../${mask}/data/new_evolution_models"
if [[ ! -e ${curr_dir5} ]]
then
        mkdir ${curr_dir5};
fi

curr_dir6="../${mask}/data/new_raxmlng_trees"
if [[ ! -e ${curr_dir6} ]]
then
        mkdir ${curr_dir6};
fi



# Creating phylogenetic trees
if [ "$radiobutton" = "prank" ]
then
        for f in $(ls ../${mask}/data/evolution_models/*prank_modeltest.out); do
                value=$(cat ${f} | grep "raxml-ng" | head -1 | awk 'ORS=" " { print $6 }');
                g=${f#../${mask}/data/evolution_models/};
                echo "../${mask}/data/multialignments/${g%_modeltest.out}.best.fas"
                raxml-ng --msa ../${mask}/data/multialignments/${g%_modeltest.out}.best.fas --model ${value} --threads $threads --prefix ../${mask}/data/raxmlng_trees/${g%_prank_modeltest.out};
        done
elif [ "$radiobutton" = "muscle" ]
then
        for f in $(ls ../${mask}/data/evolution_models/*muscle_modeltest.out); do
                value=$(cat ${f} | grep "raxml-ng" | head -1 | awk 'ORS=" " { print $6 }');
                g=${f#../${mask}/data/evolution_models/};
                raxml-ng --msa ../${mask}/data/multialignments/${g%_muscle_modeltest.out}.afa --model ${value} --threads $threads --prefix ../${mask}/data/raxmlng_trees/${g%_muscle_modeltest.out};
        done
fi

