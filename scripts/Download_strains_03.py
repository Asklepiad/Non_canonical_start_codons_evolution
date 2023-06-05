#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import json
import argparse
import shutil
import os
import sys
import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from tqdm import tqdm
from collections import defaultdict


parser = argparse.ArgumentParser()

parser.add_argument("file", type=argparse.FileType("r"), nargs="*", default=sys.stdin)
parser.add_argument("mail", type=str)

arguments = parser.parse_args()

all_species = arguments.file
Entrez.email = arguments.mail


with open (all_species, "r") as all_species_file:
    species = all_species_file.read()

all_species = species.split(", ")


# Count number of strains for each bacteria

count_complete_genome = {}
count_scaffolds = {}
id_dict = {}

for spec in tqdm(all_species):
    search_handle = Entrez.esearch(db="assembly", term=spec)
    search_record = Entrez.read(search_handle)
    count = int(search_record["Count"])

    search_handle_complete = Entrez.esearch(db="assembly", term=f'{spec}[orgn] AND ("complete genome"[filter] AND "latest refseq"[filter])', retmax=count)
    search_record_complete = Entrez.read(search_handle_complete)
    count_complete_genome[spec] = len(search_record_complete["IdList"])

    search_handle_scaffold = Entrez.esearch(db="assembly", term=f'{spec}[orgn] AND ("scaffold level"[filter] AND "latest refseq"[filter])', retmax=count)
    search_record_scaffold = Entrez.read(search_handle_scaffold)
    count_scaffolds[spec] = len(search_record_scaffold["IdList"])

    search_handle_tax = Entrez.esearch(db="taxonomy", term=spec)
    search_record_tax = Entrez.read(search_handle_tax)
    if len(search_record_tax['IdList']) == 1:
        id_dict[spec] = search_record_tax['IdList'][0]
    else:
        id_dict[spec] = 'NA'



# Create dataframe with number of strains

dd = defaultdict(list)
for all_dicts in (id_dict, count_complete_genome, count_scaffolds):
    for key, value in all_dicts.items():
        dd[key].append(value)

assembly_df = pd.DataFrame.from_dict(dd, orient='index')
assembly_df = assembly_df.reset_index() 
assembly_df.rename(columns = {'index': 'Species',0: "Taxonomy_ID", 1: "Complete_genome", 2: "Scaffolds"}, 
                   inplace = True)

assembly_df['Complete_scaffolds'] = assembly_df['Complete_genome'] + assembly_df['Scaffolds']



# Divide all bacteria into 3 groups

many_complete_genomes = [] # the bacterium has more than 100 complete genome assemblies
complete_scaffolds_100 = [] # the bacterium has less than 100 complete genome assemblies, but more than 100 complete genome plus scaffold assemblies
complete_scaffolds_less = [] # the bacterium has less than 100 complete genome plus scaffold assemblies

many_complete_genomes_id = [] 
complete_scaffolds_100_id = []
complete_scaffolds_less_id = []

assembly_df.apply(lambda x: many_complete_genomes_id.append(x['Taxonomy_ID']) if (x['Complete_genome'] > 100) else None, axis=1)
assembly_df.apply(lambda x: many_complete_genomes.append(x['Species']) if (x['Complete_genome'] > 100) else None, axis=1)
assembly_df.apply(lambda x: complete_scaffolds_100_id.append(x['Taxonomy_ID']) if ((x['Complete_genome'] < 100) & (x['Complete_scaffolds'] > 100)) else None, axis=1)
assembly_df.apply(lambda x: complete_scaffolds_100.append(x['Species']) if ((x['Complete_genome'] < 100) & (x['Complete_scaffolds'] > 100)) else None, axis=1)
assembly_df.apply(lambda x: complete_scaffolds_less_id.append(x['Taxonomy_ID']) if ((x['Complete_scaffolds'] > 20) & (x['Complete_scaffolds'] < 100)) else None, axis=1)
assembly_df.apply(lambda x: complete_scaffolds_less.append(x['Species']) if ((x['Complete_scaffolds'] > 20) & (x['Complete_scaffolds'] < 100)) else None, axis=1)



# Download genomes

# Dict of ID lists by each bacteria:
id_lists_dict = {}

# Complete genome + scaffolds < 100
for bact in complete_scaffolds_less:
    search_handle_small = Entrez.esearch(db="assembly", term= bact)
    search_record_small = Entrez.read(search_handle_small)
    count_small = int(search_record_small["Count"])

    search_handle_small_2 = Entrez.esearch(db="assembly", term= f'{bact}[orgn] AND (("complete genome"[filter] 
                                         OR "scaffold level"[filter]) AND "latest refseq"[filter])', retmax=count_small)
    search_record_small_2 = Entrez.read(search_handle_small_2)
    id_lists_dict[bact] = search_record_small_2['IdList']



# Complete genome + scaffolds > 100
for bact_2 in complete_scaffolds_100:
    search_handle_middle = Entrez.esearch(db="assembly", term= bact_2)
    search_record_middle = Entrez.read(search_handle_middle)
    count_middle = int(search_record_middle["Count"])

    # Complete genomes
    search_handle_middle_1 = Entrez.esearch(db="assembly", term= f'{bact_2}[orgn] AND ("complete genome"[filter] 
                                   AND "latest refseq"[filter])', retmax=count_middle)
    search_record_middle_1 = Entrez.read(search_handle_middle_1)
    id_lists_dict[bact_2] = search_record_middle_1['IdList']

  # Scaffolds
    count_scaff = 100 - assembly_df[assembly_df['Species'] == bact_2]['Complete_genome']
    search_handle_middle_2 = Entrez.esearch(db="assembly", term= f'{bact_2}[orgn] AND ("scaffold level"[filter] 
                                   AND "latest refseq"[filter])', retmax=count_scaff)
    search_record_middle_2 = Entrez.read(search_handle_middle_2)
    id_lists_dict[bact_2].append(search_record_middle_2['IdList'])

print(id_lists_dict)



# Bacteria with more than 100 complete genome assemblies were filtered out by the PanACoTA pipeline (to remove close assemblies by Mash genetic distance) and 100 random assemblies were taken

# Read json file with 100 random strains for each bacteria
with open("file_complete_genomes.json", "r") as file_complete:
    file_results_cut = json.load(file_complete)



# Dict: bacteria - list of esearch results for each bacteria
assembly_100_dict = {}

for key, value in tqdm(file_results_cut.items()):
    bact_id = key.split('-')[1]
    bact_name = '_'.join(assembly_df[assembly_df['Taxonomy_ID'] == bact_id]['Species'].values[0].split(' '))
    assembly_100_dict[bact_name] = []

    for strain in value:
        strain_name = strain.split('tmp_files/')[1].split('.')[0]

        search_handle_1 = Entrez.esearch(db="assembly", term= strain_name)
        search_record_1 = Entrez.read(search_handle_1)

        assembly_100_dict[back_name].append(search_record_1)



# Add all bacteria to dict:

for key, value in assembly_100_dict.items():
    id_lists_dict[key] = []

    for res in value:
        if len(res['IdList']) == 1:
            id_lists_dict[key].append(res['IdList'][0])
        else:
            id_lists_dict[key].append(res['IdList'][0])
            id_lists_dict[key].append(res['IdList'][1])



# Save dict in json format
with open("id_lists.json", "w") as f:
    json.dump(id_lists_dict, f)

print('All strains were downloaded and json file with ID lists was saved in current directory!')

