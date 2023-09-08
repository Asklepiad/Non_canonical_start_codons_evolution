#!/usr/bin/env python
# coding: utf-8



import Bio
from Bio import Entrez
from tqdm import tqdm

import json
import os
import argparse
import sys
import random


# Arguments

parser = argparse.ArgumentParser()

parser.add_argument("file", type=str, default=sys.stdin) # argparse.FileType("r")
parser.add_argument("mail", type=str)


arguments = parser.parse_args()

all_species = arguments.file
Entrez.email = arguments.mail



# Read file

with open(all_species, "r") as f:
    species = json.load(f)

current_species = []

for key, value in species.items():
    current_species.append(key)

print('All species:')
print(current_species)



# Count number of complete genome, scaffold and chromosome assemblies

count_complete_genome = {}
count_scaffolds = {}
count_chromosome = {}

ids_complete_genome = {}
ids_scaffolds = {}
ids_chromosome = {}

count_complete_scaffold_chromosome = {}



for spec in tqdm(current_species):
    search_handle_1 = Entrez.esearch(db="assembly", term= spec)
    search_record_1 = Entrez.read(search_handle_1)
    count = int(search_record_1["Count"])

    # Complete genome
    search_handle_2 = Entrez.esearch(db="assembly", term= f'{spec}[Organism] AND ("complete genome"[Assembly Level])', retmax=count)
    search_record_2 = Entrez.read(search_handle_2)
    count_complete_genome[spec] = len(search_record_2["IdList"])
    ids_complete_genome[spec] = search_record_2["IdList"]
    
    # Scaffolds
    search_handle_3 = Entrez.esearch(db="assembly", term= f'{spec}[Organism] AND ("scaffold"[Assembly Level])', retmax=count)
    search_record_3 = Entrez.read(search_handle_3)
    count_scaffolds[spec] = len(search_record_3["IdList"])
    ids_scaffolds[spec] = search_record_3["IdList"]

    # Chromosome
    search_handle_4 = Entrez.esearch(db="assembly", term= f'{spec}[Organism] AND ("chromosome"[Assembly Level])', retmax=count)
    search_record_4 = Entrez.read(search_handle_4)
    count_chromosome[spec] = len(search_record_4["IdList"])
    ids_chromosome[spec] = search_record_4["IdList"]
    
    count_complete_scaffold_chromosome[spec] = len(search_record_2["IdList"]) + len(search_record_3["IdList"])+ len(search_record_4["IdList"])



# Divide all bacteria into 2 groups

bact_for_panacota = []
bact_less_50_complete = []

for key, value in count_complete_genome.items():
    if value > 50:
        bact_for_panacota.append(key)
    else:
        bact_less_50_complete.append(key)


# Divide second group of bacteria into 2 groups
panacota_result_dict = {}

for bact in bact_less_50_complete:
    
    if count_complete_scaffold_chromosome[bact] >= 50:
        
        if count_complete_genome[bact] + count_scaffolds[bact] >= 50:
            count_scaff = 50 - count_complete_genome[bact]
            panacota_result_dict[bact] = ids_complete_genome[bact] + random.sample(ids_scaffolds[bact], k=count_scaff, )
    
        else:
            count_chrom = 50 - count_complete_genome[bact] - count_scaffolds[bact]
            panacota_result_dict[bact] = ids_complete_genome[bact] + ids_scaffolds[bact] + random.sample(ids_chromosome[bact], k=count_chrom, )

    else:
        panacota_result_dict[bact] = ids_complete_genome[bact] + ids_scaffolds[bact] + ids_chromosome[bact]





print('Bacteria for panacota')
print(bact_for_panacota)




# Run Panacota pipeline
os.system("mkdir Results_panacota")
os.chdir("Results_panacota")

for big_bact in bact_for_panacota:
    os.system(f'PanACoTA prepare -g "{big_bact}" -l complete')





# Dictionary with assemblies after panacota
file_results = {}
num_bact = -1

for root, dirs, files in os.walk("./"):
    for file in files:
        if file.startswith('assembly'):
            num_bact += 1
            # print(file)
            # print(num_bact)
            with open("./" + file.split('-')[1].split('.')[0] + "/" + file, 'r') as result_file:
                file_results[file] = result_file.readlines()
                




# Delete first row with heading
for key, value in file_results.items():
    del value[0]



# Make a dictionary with 50 random assemblies

random.seed(10)
file_results_cut = {}

for key, value in file_results.items():
    if len(value) > 50:
        file_results_cut[key] = random.sample(value, k=50, )
    else:
        file_results_cut[key] = value




for key, value in file_results_cut.items():
    panacota_result_dict[key.split('-')[1].split('.')[0]] = []
    for assembl in value:
        panacota_result_dict[key.split('-')[1].split('.')[0]].append(assembl.split('.')[0])


# Save a dictionary with cut number of assemblies

with open("file_filtered_genomes.json", "w") as f:
    json.dump(panacota_result_dict, f)






