
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

arguments = parser.parse_args()
path_to_alignments = arguments.path

list_align_files = os.listdir(path_to_alignments)

all_genes_in_strain = {}
strain_gap_dict = {}
bad_strains_dict = {}

for align_file in list_align_files:
  alignment = AlignIO.read(open(path_to_alignments + align_file), "fasta")
  alignment_length = alignment.get_alignment_length()

  atg_cnt = 0
  gtg_cnt = 0
  ttg_cnt = 0
  gap_cnt = 0

  for record in alignment: 

    # Write number of genes in each strain
    if record.id[15:] in all_genes_in_strain:
      all_genes_in_strain[record.id[15:]] += 1
    else:
      all_genes_in_strain[record.id[15:]] = 1

    # Check gaps
    if record.seq[:3] == '---':
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
for d in (strain_gap_dict, all_genes_in_strain, bad_strains_dict):
    for key, value in d.items():
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
  strain_df.iloc[row, 4] = strain_df.iloc[row, 1] / strain_df.iloc[row, 2] * 100
  strain_df.iloc[row, 5] = strain_df.iloc[row, 3] / strain_df.iloc[row, 2] * 100

# Medians and std
med_gaps = strain_df['Proportion_gaps'].median()
med_minor = strain_df['Proportion_minor'].median()
std_gaps = np.std(strain_df['Proportion_gaps'])
std_minor = np.std(strain_df['Proportion_minor'])

for row in range(len(strain_df)):
  if strain_df.iloc[row, 4] >= (med_gaps + 2 * std_gaps):
    list_many_gaps_strain.append(strain_df.iloc[row, 0])
  if strain_df.iloc[row, 5] >= (med_minor + 2 * std_minor):
    list_many_minors_strain.append(strain_df.iloc[row, 0])

# Barplots
sns.set(rc={'figure.figsize':(13,6)})
fig, axs = plt.subplots(ncols=2)
strain_plot_1 = sns.barplot(data=strain_df, x = 'Strain', y="Proportion_gaps", 
                        order=strain_df.sort_values('Proportion_gaps').Strain, ax=axs[0])
strain_plot_1.set_xticklabels(strain_plot_1.get_xticklabels(), rotation=90);


strain_plot_2 = sns.barplot(data=strain_df, x = 'Strain', y="Proportion_minor", 
                        order=strain_df.sort_values('Proportion_minor').Strain, ax=axs[1])
strain_plot_2.set_xticklabels(strain_plot_2.get_xticklabels(), rotation=90);  

return list_many_gaps_strain, list_many_minors_strain
