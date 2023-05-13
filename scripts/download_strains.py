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


#Entrez.email = "valeriyakudr98@gmail.com"
Entrez.api_key = "a676bc7d097855d36e502959d6ef43aa4e08"



#generalists = ['Acinetobacter baumannii', 'Acinetobacter indicus', 'Agrobacterium tumefaciens', 'Bacillus cereus',
#               'Bacillus megaterium', 'Bacillus pumilus', 'Bacillus subtilis', 'Bacillus thuringiensis', 
#               'Clostridium botulinum', 'Enterococcus faecalis', 'Escherichia coli', 'Klebsiella pneumoniae', 
#               'Lactococcus lactis', 'Mycobacterium avium', 'Pseudomonas aeruginosa', 'Pseudomonas fluorescens', 
#               'Pseudomonas putida', 'Stenotrophomonas maltophilia']

#specialists = ['Bacillus altitudinis', 'Bifidobacterium bifidum', 'Burkholderia cenocepacia', 
#               'Candidatus Carsonella ruddii', 'Chlamydia trachomatis', 'Corynebacterium diphtheriae', 
#               'Corynebacterium pseudotuberculosis', 'Cupriavidus taiwanensis', 'Enterobacter hormaechei',
#               'Francisella tularensis', 'Histophilus somni', 'Levilactobacillus brevis', 'Lactobacillus johnsonii', 
#               'Lactobacillus salivarius', 'Leptospira interrogans', 'Mesoplasma florum', 'Morganella morganii', 
#               'Mycoplasma bovis', 'Neisseria gonorrhoeae', 'Neisseria meningitidis', 'Paenibacillus larvae', 
#               'Pantoea ananatis', 'Pediococcus acidilactici', 'Porphyromonas gingivalis', 
#               'Propionibacterium freudenreichii', 'Rhodopseudomonas palustris', 'Riemerella anatipestifer',
#               'Salinibacter ruber', 'Staphylococcus haemolyticus', 'Staphylococcus simulans', 'Vibrio anguillarum',
#               'Vibrio campbellii', 'Vibrio cholerae', 'Weissella cibaria', 'Xanthomonas campestris', 
#               'Xylella fastidiosa', 'Zymomonas mobilis']

#all_species = generalists + specialists

parser = argparse.ArgumentParser()

parser.add_argument("file", type=argparse.FileType("r"), nargs="*", default=sys.stdin)

arguments = parser.parse_args()

all_species = []
for texts in arguments.file:
        text = texts.read()
        result_list = text.strip().split(", ")
        all_species += (result_list)



# Count number of strains for each bacteria

count_complete_genome = {}
count_scaffolds = {}
id_dict = {}

for spec in tqdm(all_species):
  search_handle_1 = Entrez.esearch(db="assembly", term=spec)
  search_record_1 = Entrez.read(search_handle_1)
  count = int(search_record_1["Count"])

  search_handle_2 = Entrez.esearch(db="assembly", term=f'{spec}[orgn] AND ("complete genome"[filter] AND "latest refseq"[filter])', retmax=count)
  search_record_2 = Entrez.read(search_handle_2)
  count_complete_genome[spec] = len(search_record_2["IdList"])

  search_handle_2 = Entrez.esearch(db="assembly", term=f'{spec}[orgn] AND ("scaffold level"[filter] AND "latest refseq"[filter])', retmax=count)
  search_record_2 = Entrez.read(search_handle_2)
  count_scaffolds[spec] = len(search_record_2["IdList"])

  search_handle_3 = Entrez.esearch(db="taxonomy", term=spec)
  search_record_3 = Entrez.read(search_handle_3)
  if len(search_record_3['IdList']) == 1:
    id_dict[spec] = search_record_3['IdList'][0]
  else:
    id_dict[spec] = 'NA'


# Create dataframe with number of strains

dd = defaultdict(list)
for d in (id_dict, count_complete_genome, count_scaffolds):
    for key, value in d.items():
        dd[key].append(value)

assembly_df = pd.DataFrame.from_dict(dd, orient='index')
assembly_df = assembly_df.reset_index() 
assembly_df.rename(columns = {'index': 'Species',0: "Taxonomy_ID", 1: "Complete_genome", 2: "Scaffolds"}, inplace = True)

assembly_df['Complete_scaffolds'] = 0
for row in range(len(assembly_df)):
  assembly_df.iloc[row, 4] = assembly_df.iloc[row, 2] + assembly_df.iloc[row, 3]


# Divide all bacteria into 3 groups

many_complete_genomes_id = []
complete_scaffolds_100_id = []
complete_scaffolds_less_id = []

many_complete_genomes = []
complete_scaffolds_100 = []
complete_scaffolds_less = []

for row in range(len(assembly_df)):
  if assembly_df.iloc[row, 2] > 100:
    many_complete_genomes_id.append(assembly_df.iloc[row, 1])
    many_complete_genomes.append(assembly_df.iloc[row, 0])
  else:
    if assembly_df.iloc[row, 4] > 100:
      complete_scaffolds_100_id.append(assembly_df.iloc[row, 1])
      complete_scaffolds_100.append(assembly_df.iloc[row, 0])
    else:
      if assembly_df.iloc[row, 4] > 20:
        complete_scaffolds_less_id.append(assembly_df.iloc[row, 1])
        complete_scaffolds_less.append(assembly_df.iloc[row, 0])


# Download genomes

# Dict of ID lists by each bacteria:

id_lists_dict = {}

# Complete genome + scaffolds < 100
for bact in complete_scaffolds_less:
  search_handle_1 = Entrez.esearch(db="assembly", term= bact)
  search_record_1 = Entrez.read(search_handle_1)
  count = int(search_record_1["Count"])

  search_handle_2 = Entrez.esearch(db="assembly", term= f'{bact}[orgn] AND (("complete genome"[filter] OR "scaffold level"[filter]) AND "latest refseq"[filter])', retmax=count)
  search_record_2 = Entrez.read(search_handle_2)
  id_lists_dict[bact] = search_record_2['IdList']


# Complete genome + scaffolds > 100
for bact_2 in complete_scaffolds_100:
  search_handle_1 = Entrez.esearch(db="assembly", term= bact_2)
  search_record_1 = Entrez.read(search_handle_1)
  count = int(search_record_1["Count"])

  # Complete genomes
  search_handle_3 = Entrez.esearch(db="assembly", term= f'{bact_2}[orgn] AND ("complete genome"[filter] AND "latest refseq"[filter])', retmax=count)
  search_record_3 = Entrez.read(search_handle_3)
  id_lists_dict[bact_2] = search_record_3['IdList']

  # Scaffolds
  count_scaff = 100 - assembly_df[assembly_df['Species'] == bact_2]['Complete_genome']
  search_handle_4 = Entrez.esearch(db="assembly", term= f'{bact_2}[orgn] AND ("scaffold level"[filter] AND "latest refseq"[filter])', retmax=count_scaff)
  search_record_4 = Entrez.read(search_handle_4)
  id_lists_dict[bact_2].append(search_record_4['IdList'])

print(id_lists_dict)
# After panacota

# Read json file with 100 random strains for each bacteria
with open("file_complete_genomes.json", "r") as f:
  file_results_cut = json.load(f)

# Dict: bacteria - list of esearch results for each bacteria
assembly_100_dict = {}

for key, value in tqdm(file_results_cut.items()):
  bact_id = key.split('-')[1]
  back_name = '_'.join(assembly_df[assembly_df['Taxonomy_ID'] == bact_id]['Species'].values[0].split(' '))
  assembly_100_dict[back_name] = []
  # print(back_name)
  # print(len(value))

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
