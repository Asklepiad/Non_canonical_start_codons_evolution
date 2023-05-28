#!/bin/bash


# Variables
path=$1
aligner=$2
po_parameter=$3
mail=$4
threads=$5
pre_short=${path#../data/jsons/}
short_name=${pre_short%.json}
echo "$short_name" "bu-ga-ga"

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
conda run -n python_start_codons python3 Parsing_NCBI_1.py "$path" "$mail"
echo "script 1 done"
conda run -n prokka_start_codons ./prokka_annotate2.sh "$short_name"
echo "script 2 done"
conda run -n python_start_codons python3 First_table_creating3.py "$short_name"
echo "script 3 done"
conda run -n proteinortho_start_codons ./proteinortho_script4.sh "$short_name" "$po_parameter"
echo "script 4 done"
conda run -n python_start_codons python3 Muscle_preparing_5.py "$path"
echo "script 5 done"
conda run -n R_start_codons Rscript Statscript.R "$short_name"
echo "statistics script done"
if [ "$aligner" = "prank" ]
then
	conda run -n prank_start_codons ./prank_align6.sh "$short_name"
elif [ "$aligner" = "muscle" ]
then
	conda run -n muscle_start_codons ./muscle_align6.sh "$short_name" "$threads"
fi
echo "script 6 done"
conda run -n modeltest_start_codons ./modeltree_maker7.sh "$short_name" "$aligner"
echo "script 7 done"
conda run -n raxml_start_codons ./raxml_tree8.sh "$short_name" "$aligner" "$threads"
echo "script 8 done"
conda run -n python_start_codons python3 Config_tree_maker_9.py "$short_name" "$aligner"
echo "script 9 done"
conda run -n python_start_codons python3 Ete3_maker_10.py "$short_name"
echo "script 10 done"
