# !/bin/bash

#mkdir ../raxmlng_trees
mask=${C_psittaci}
for f in $(ls ${mask}*.out); do
	value=$(cat ${f} | grep "raxml-ng" | head -1 | awk 'ORS=" " { print $6 }');
	raxml-ng --msa ../multialignments/${f%modeltest.out}.afa --model ${value} --threads 7 --prefix ../raxmlng_trees/${f%modeltest.out};
done
