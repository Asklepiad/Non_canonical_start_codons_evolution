#!/bin/bash


# Reads arguments
input=$1
radiobutton=$2


# Choosing the best evolution model for phylogenetic tree
if [ "$radiobutton" = "prank" ]
then
        for f in $(ls ../${input}/data/multialignments/*_prank.best.fas); do
                echo "prank"
                echo $f
                g=${f#../${input}/data/multialignments/};
                modeltest-ng static -i ${f} -o ../${input}/data/evolution_models/${g%.best.fas}_modeltest;
        done
elif [ "$radiobutton" = "muscle" ]
then
        for f in $(ls ../${input}/data/multialignments/*.afa); do
                echo "muscle"
                echo $f
                g=${f#../${input}/data/multialignments/};
                modeltest-ng static -i ${f} -o ../${input}/data/evolution_models/${g%.afa}_muscle_modeltest;
        done
fi

