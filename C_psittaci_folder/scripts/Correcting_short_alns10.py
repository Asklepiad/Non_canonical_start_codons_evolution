#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import shutil
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from itertools import chain
from tqdm import tqdm
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

parser = argparse.ArgumentParser()

parser.add_argument("folder_name", type=str)
parser.add_argument("aligner", type=str)

arguments = parser.parse_args()

folder_name = arguments.folder_name
aligner = arguments.aligner


if aligner == "prank":
    last_3 = "fas"
elif aligner == "muscle":
    last_3 = "afa"


current_path = f"../{folder_name}/data/multialignments/"
alignments = os.listdir(current_path)       # Create list of directory content
alignments = [alignment for alignment in alignments if alignment[-3:] == last_3]
first_table = pd.read_csv(f"../{folder_name}/data/First_table.csv")   # Load first table with data
fasta_path = f"../{folder_name}/data/for_prokka_fasta/"       
hand_assemblies = os.listdir(fasta_path)     # Assemblies for "long" alignments


pattern = r"^[-]*"
orto_pattern = r"_[\d]+[\D]"
orto_rows, new_sc, new_als, genes = [], [], [], []
for alignment in tqdm(alignments):   # Looping all non-uniform alignments
    number = re.search(orto_pattern, alignment)
    orto_row = int(number.group()[1:-1])
    cur_path = os.path.join(current_path, alignment)
    fasta_iter = SeqIO.parse(cur_path, "fasta")
    
    for fasta in fasta_iter:       
      
        if fasta.seq[:3] == "---":
            string = fasta.seq
            length = re.match(pattern, str(string)).end()
            name, source = fasta.name[:14], fasta.name[15:-4]
            start = first_table.query("id == @name")["first"].iloc[0]
            finish = first_table.query("id == @name")["last"].iloc[0]
            strand = first_table.query("id == @name")["strand"].iloc[0]
            new_start = start - length
            current_source = os.path.join(fasta_path, f"{source}.fasta")
            sourcer = SeqIO.read(current_source, "fasta").seq
                        
            if new_start >= 0:
                if strand == -1:
                    new_al = sourcer.complement()[new_start:finish][::-1]
                elif strand == 1:
                    new_al = sourcer[new_start:finish]
            else:
                if strand == -1:
                    new_al = (sourcer.complement()[new_start:] + sourcer.complement()[:finish])[::-1]
                elif strand == 1:
                    new_al = sourcer[new_start:] + sourcer[:finish]
                    
            new_als.append(str(new_al))
            new_sc.append(str(new_al[:3]))
            genes.append(name)
            
# Translating new aa_seqs
new_aas = []
pat = r"^[A-Z]*\*"   # Trimming after first stop codone
for seq in new_als:
    if seq[0:3] in ["ATG", "GTG", "TTG"]:
        new_aas.append(re.search(pat, str("M" + Seq(seq[3:]).translate())).group()[:-1])
    else:
        new_aas.append("")
new_aas

# Creating table with new starts
alt_nonunif_sc = pd.DataFrame(zip(genes, new_sc, new_als, new_aas))
alt_nonunif_sc.rename(columns = {
    0: "id",
    1: "new_sc",
    2: "new_als",
    3: "new_aas"
}, inplace=True)

summary_rows = pd.read_csv(f"../{folder_name}/data/summary_rows_prokka.csv")

# Creating cog dictionary
cog_dict = {}
current_path = "../data/bigcog/"
cog_cats = os.listdir(current_path)
for cog_cat in tqdm(cog_cats):
    table = pd.read_csv(os.path.join(current_path, cog_cat), sep="\t")
    if len(table.Cat):
        cat = table.Cat.tolist()
        cat_ful = list(map(lambda x: x.split(" "),cat))
        cat_ful_single = list(chain.from_iterable(cat_ful))
        key = max(cat_ful_single, key = cat_ful_single.count)
        value = table.COG.tolist()
        cog_dict[key] = value
cog_dict["S"].append("absent")
dfft = pd.read_csv(f"../{folder_name}/data/First_table.csv")
cog = dfft["cog"].to_list()
for key in tqdm(cog_dict.keys()):   # Passing over COG categories
    cog_list = []             # For each category create list
    for cog_id in cog:        # Passing out every item
        if cog_id not in cog_dict[key]:     # For every gene, if gene's id not in dict
            cog_list.append(0)              # Mark it
        else:                               
            cog_list.append(1)              # Else mark another
    summary_rows[key] = pd.Series(cog_list)   # Creating 25 new columns

summary_rows.drop(["new_sc", "new_als", "new_aas"], axis=1, inplace=True)
alt_sc = summary_rows.merge(alt_nonunif_sc, on="id", how="outer")
alt_sc["new_sc"].fillna(value=alt_sc["start_codone"], inplace=True)
alt_sc["new_als"].fillna(value=alt_sc["n_sequence"], inplace=True)
alt_sc["new_aas"].fillna(value=alt_sc["aa_sequence"], inplace=True)
alt_sc["is_new"] = alt_sc["n_sequence"] == alt_sc["new_als"]
alt_sc["new_sc_c"] = [i[:3] for i in alt_sc["new_als"]]

yes_no = []   # Is n_sequence has equal length with a_sequence?
for i in range(len(alt_sc)):
    delta = len(alt_sc["new_als"][i]) - (len(alt_sc["new_aas"][i]) * 3 + 3)
    if delta > 0:
        alt_sc["new_als"][i] = alt_sc["new_als"][i][:-delta]
    if len(alt_sc["new_als"][i]) == 3:
        alt_sc["new_als"][i] = alt_sc["n_sequence"][i]
        alt_sc["new_aas"][i] = alt_sc["aa_sequence"][i]
    yes_no.append(len(alt_sc["new_als"][i]) - (len(alt_sc["new_aas"][i]) * 3 + 3))
alt_sc.to_csv(f"../{folder_name}/data/summary_rows_new.csv", index=False)


# Now we want to recompute multialignments.

non_uniform_or_list = list(set(alt_sc.query("uniformity == 'non-uniform'").ortologus_row))
for orto_row in tqdm(non_uniform_or_list):
    subset = alt_sc.query("ortologus_row == @orto_row")
    with open (f"../{folder_name}/data/new_multialignments/{folder_name}_withoutp_{str(orto_row)}.fasta", "w") as nucleotide_fasta:
        for index, row in subset.iterrows():
            if row['type_of_the_gene'] != "pseudogene":
                nucleotide_fasta.write(">")
                nucleotide_fasta.write(row["id"])
                nucleotide_fasta.write("_")
                nucleotide_fasta.write(row["source"])
                nucleotide_fasta.write("_")
                nucleotide_fasta.write(row["new_sc_c"])
                nucleotide_fasta.write("\n")
                nucleotide_fasta.write(row["new_als"])
                nucleotide_fasta.write("\n")

