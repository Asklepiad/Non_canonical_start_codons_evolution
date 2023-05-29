#!/bin/bash


# Reads arguments
proj=$1
identity=$2


# Creates folder for storing results of multialignments
curr_dir="../${proj}/data/multialignments"
if [[ ! -e ${curr_dir} ]]
then
        mkdir ${curr_dir};
fi


# Construction of ortologous groups
proteinortho6.pl --project=${proj} ../${proj}/data/orto_rows/*.fasta --debug --selfblast --singles --identity=${identity}
mv ${proj}.proteinortho* ${proj}.info ${proj}.blast* ../${proj}/data/
