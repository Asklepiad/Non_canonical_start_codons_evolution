#!/bin/bash


# Variables
organism=$1
aligner=$2
po_parameter=$3
mail=$4
first_name=${organism:0:1}"_"
presecond_name=$(echo "$organism" | awk -F " " '{print $2}')
second_name=${presecond_name:0:9}
short_name=${first_name}${second_name}


# Folders
input="../${short_name}"
if [[ ! -e ${input} ]]
then
	mkdir ${input};
	mkdir "${input}/data"
	mkdir "${input}/figures"
	mkdir "${input}/figures/trees"
fi

curr_dir="${input}/data/for_prokka_fasta"
if [[ ! -e ${curr_dir} ]]
then
        mkdir ${curr_dir};
fi


# Scripts
echo "script 0 done"
conda run -n base python3 Parsing_NCBI_1.py "$organism" "$mail"
echo "script 1 done"
conda run -n prokka_env ./prokka_annotate2.sh "$short_name"
echo "script 2 done"
conda run -n base python3 First_table_creating3.py "$short_name"
echo "script 3 done"
conda run -n proteinortho ./proteinortho_script4.sh "$short_name" "$po_parameter"
echo "script 4 done"
conda run -n base python3 Muscle_preparing_5.py "$organism"
echo "script 5 done"
conda run -n base Rscript Statscript.R "$short_name"
echo "statistics script done"
conda run -n multialign ./prank_align6.sh "$short_name"
echo "script 6 done"
conda run -n modeltest ./modeltree_maker7.sh "$short_name" "$aligner"
echo "script 7 done"
conda run -n raxml ./raxml_tree8.sh "$short_name" "$aligner"
echo "script 8 done"
conda run -n base python3 Config_tree_maker_9.py "$short_name" "$aligner"
echo "script 9 done"
conda run -n base python3 Ete3_maker_10.py "$short_name"
echo "script 10 done"
#conda run -n base python3 Correcting_short_alns10.py "$short_name" "$aligner"
#echo "script 10 done"
