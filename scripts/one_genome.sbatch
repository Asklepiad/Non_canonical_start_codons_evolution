#!/bin/bash

#SBATCH --job-name=s1ta1rt_codons
#SBATCH --error=q2new_start_codons.err
#SBATCH --output=q2new_start_codons.log
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --partition=all

time(./folder_creators_server0.sh "S_ruber.json" "muscle" 75 "bogdan.sotnikov.1999@mail.ru" 24)
