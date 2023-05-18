#!/bin/bash

#SBATCH --job-name=start_codons
#SBATCH --error=start_codons.err
#SBATCH --output=start_codons.log
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --partition=all

#conda run -n python_start_codons python3 Ete3_maker_10.py "F_varium"
#conda run -n R_start_codons Rscript Statscript.R "F_varium"
time(./json.sh "../data/id_lists.json" "muscle" 75 "bogdan.sotnikov.1999@mail.ru" 24)
