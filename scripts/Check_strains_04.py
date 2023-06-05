#!/usr/bin/env python
# coding: utf-8


from collections import defaultdict
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("path", type=str)
parser.add_argument("aligner", type=str)

arguments = parser.parse_args()
path_to_alignments = arguments.path
aligner = arguments.aligner


if aligner == "prank":
    num1 = 9
    word = ".best.fas"
elif aligner == "muscle":
    num1 = 4
    word = ".afa"
list_align_files = [align for align in os.listdir(path_to_alignments) if align[-num1:] == word]


all_genes_in_strain = {}
strain_gap_dict = {}
bad_strains_dict = {}


for align_file in list_align_files:
    alignment = AlignIO.read(f"{path_to_alignments}{align_file}", "fasta")
    alignment_length = alignment.get_alignment_length()

    atg_cnt = 0
    gtg_cnt = 0
    ttg_cnt = 0
    gap_cnt = 0

    for record in alignment: 

    # Write number of genes in each strain
        if record.id[15:] in all_genes_in_strain: # the first 15 letters in the name of the gene is the identifier of the gene, after it comes the name of the strain
            all_genes_in_strain[record.id[15:]] += 1
        else:
            all_genes_in_strain[record.id[15:]] = 1

    # Check gaps
        if record.seq[:3] == '---': # the first 3 letters in the alignment is the start codon
            if record.id[15:] in strain_gap_dict:
                strain_gap_dict[record.id[15:]] += 1
            else:
                strain_gap_dict[record.id[15:]] = 1 


    # Check minor start codons
        if record.seq[:3] == 'ATG':
            atg_cnt += 1
        elif record.seq[:3] == 'GTG':
            gtg_cnt += 1
        elif record.seq[:3] == 'TTG':
            ttg_cnt += 1
        elif record.seq[:3] == '---':
            gap_cnt += 1

    var = {atg_cnt:"atg_cnt", gtg_cnt:"gtg_cnt", ttg_cnt:"ttg_cnt", gap_cnt:"gap_cnt"}
    major = var.get(max(var))[:3].upper()


    for record in alignment:
        if record.seq[:3] != major:
            if record.id[15:] in bad_strains_dict:
                bad_strains_dict[record.id[15:]] += 1
            else:
                bad_strains_dict[record.id[15:]] = 1



# Combining dicts
dd = defaultdict(list)
for all_dicts in (strain_gap_dict, all_genes_in_strain, bad_strains_dict):
    for key, value in all_dicts.items():
        dd[key].append(value)

        

# Make dataframe
strain_df = pd.DataFrame.from_dict(dd, orient='index')
strain_df = strain_df.reset_index() 

strain_df.rename(columns = {0: "Gaps", 1: "All", 'index': 'Strain', 2: "Minor_start"}, inplace = True)
strain_df['Proportion_gaps'] = 0
strain_df['Proportion_minor'] = 0 

list_many_gaps_strain = []
list_many_minors_strain = []

for row in range(len(strain_df)):
    strain_df['Proportion_gaps'] = strain_df['Gaps'] / strain_df['All'] * 100
    strain_df['Proportion_minor'] = strain_df['Minor_start'] / strain_df['All'] * 100



# Medians and std
med_gaps = strain_df['Proportion_gaps'].median()
med_minor = strain_df['Proportion_minor'].median()
std_gaps = np.std(strain_df['Proportion_gaps'])
std_minor = np.std(strain_df['Proportion_minor'])

list_many_gaps_strain = []
list_many_minors_strain = []

strain_df.apply(lambda x: list_many_gaps_strain.append(x['Strain']) if (x['Proportion_gaps'] > (med_gaps + 2 * std_gaps)) else None, axis=1)
strain_df.apply(lambda x: list_many_minors_strain.append(x['Strain']) if (x['Proportion_minor'] > (med_minor + 2 * std_minor)) else None, axis=1)



# Barplots
sns.set(rc={'figure.figsize':(13,6)})
fig, axs = plt.subplots(ncols=2)
strain_plot_1 = sns.barplot(data=strain_df, x = 'Strain', y="Proportion_gaps", 
                        order=strain_df.sort_values('Proportion_gaps').Strain, ax=axs[0])
strain_plot_1.set_xticklabels(strain_plot_1.get_xticklabels(), rotation=90);


strain_plot_2 = sns.barplot(data=strain_df, x = 'Strain', y="Proportion_minor", 
                        order=strain_df.sort_values('Proportion_minor').Strain, ax=axs[1])
strain_plot_2.set_xticklabels(strain_plot_2.get_xticklabels(), rotation=90);  

print(list_many_gaps_strain, list_many_minors_strain)

