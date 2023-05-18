#!/bin/bash

#SBATCH --job-name=start_codons
#SBATCH --error=z_start_codons.err
#SBATCH --output=z_start_codons.log
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --partition=gpu

time(conda run -n R_start_codons Rscript Statscript.R "S_ruber")
