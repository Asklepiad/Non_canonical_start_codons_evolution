#!/bin/bash


# Reads arguments
mask=$1


curr_dir1="../${mask}/data/${mask}_annotate"
if [[ ! -e ${curr_dir1} ]]
then
        mkdir ${curr_dir1};
fi


# Creates folder for storing fasta-files for further processing by proteinortho script
curr_dir2="../${mask}/data/orto_rows"
if [[ ! -e ${curr_dir2} ]]
then
        mkdir ${curr_dir2};
fi


# Reannotates each DNA molecule (chromosome or plasmid)
curr_dir3="../${mask}/data/for_prokka_fasta"
for f in $(ls ${curr_dir3}/${mask}*.fasta); do
	g=$(echo ${f#${curr_dir3}/});
        prokka --outdir ${f%.fasta} --cpus 0 --prefix ${g%.fasta} ${f};
done


# Moves annotations (gbk files) into ${mask}_annotate folder
for j in $(ls ${curr_dir3}/${mask}*); do
        if [[ "${j}" == "${j%.*}.gbk" ]]; then
		h=${j%.gbk}
                mv ${curr_dir3}/${h}/${h}.gbk ${curr_dir1}
        fi
done
