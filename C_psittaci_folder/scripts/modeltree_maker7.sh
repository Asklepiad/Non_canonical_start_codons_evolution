#!/bin/bash

input=$1
radiobutton=$2


if [[ radiobutton=="prank" ]]
then
	for f in $(ls ../${input}/data/multialignments/*_prank.best.fas); do
		g=${f#../${input}/data/multialignments/};
        	modeltest-ng static -i ${f} -o ../${input}/data/evolution_models/${g%.best.fas}_modeltest;
	done
elif [[ radiobutton=="muscle" ]]
then
	for f in $(ls ../${input}/data/multialignments/{input}*.afa); do
		g=${f#../${input}/data/multialignments/};
        	modeltest-ng static -i ${f} -o ../${input}/data/evolution_models/${g%.afa}modeltest;
	done
fi
