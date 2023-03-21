#!/usr/bin/env python
# coding: utf-8

# On the next stage I need to download file with ortologus rows in the system. I have mounted google drive with project files (in the future I'm going to use GitHub for these purposes). 

import os
import sys
import argparse
import shutil
import numpy as np
import pandas as pd
from tqdm import tqdm


parser = argparse.ArgumentParser()

parser.add_argument("organism_name", type=str)

arguments = parser.parse_args()

organism_name = arguments.organism_name
lst = organism_name.split()
folder_name = f"{lst[0][0]}_{lst[1]}"


# ### Orto rows with singletons, but without paralog-containing rows
# Uploading and modificating orto rows tables
orto_rows = pd.read_csv(f"../{folder_name}/data/{folder_name}.proteinortho.tsv", sep="\t") # Uploading dataframe with orto rows
orto_rows = orto_rows.rename(columns = {"# Species": "Strains"})
orto_rows["ortologus_row"] = orto_rows.index + 1 # Creating ortorows numbers

# Uploading first csv-table and creating a new column in it
orto_rows_list = orto_rows.index
df1 = pd.read_csv(f"../{folder_name}/data/First_table.csv")
df1["ortologus_row"] = 0

#Filling orto_row column (sounds like an oxymoron)
organism = organism_name
for strain in tqdm(set(df1["p_c_unity"])):
    maskstring = f"{folder_name}{strain}.fasta"
    orrow = orto_rows.loc[:, maskstring].str[:14]
    for index in orto_rows_list:
        df1.loc[df1[df1["id"] == orrow[index]].index, "ortologus_row"] = index + 1

# Number of the paralogs
print("Number of the paralogs =", sum(orto_rows.query("Genes > Strains").Genes) - sum(orto_rows.query("Genes > Strains").Strains))

# Creating a subset without pararows
pararows_numbers = orto_rows.query("Genes > Strains").ortologus_row
df1 = df1.query("ortologus_row not in @pararows_numbers").query("ortologus_row != 0")

# Compiling data about start-codons of ortologus rows
start_codon_per_row = df1.groupby("ortologus_row", as_index=False).agg({"start_codone": ".".join})
start_codon_per_row["start_codone"]
orr_start_list = []
for number in tqdm(range(len(start_codon_per_row))):
    orr_start_list.append(start_codon_per_row.iloc[number, 1].split("."))
start_codons = pd.DataFrame(
    {
        "ortologus_row": start_codon_per_row["ortologus_row"],
        "start_codons": pd.Series(orr_start_list),
        "ATG": 0.0,
        "GTG": 0.0,
        "TTG": 0.0,
    }
)

# Computing frequencies of exact start-codons
start_codons.sort_values("ortologus_row", inplace=True)
for row in tqdm(range(len(start_codons))):
    freqs = pd.Series(start_codons["start_codons"][row]).value_counts() # If it will be need to visualize percents, use normalize=True in brackets
    rowlength = len(start_codons["start_codons"][row])
    if "ATG" in freqs:
        start_codons["ATG"][row] = freqs["ATG"]
    else:
        start_codons["ATG"][row] = 0
    if "GTG" in freqs:
        start_codons["GTG"][row] = freqs["GTG"]
    else:
        start_codons["GTG"][row] = 0
    if "TTG" in freqs:
        start_codons["TTG"][row] = freqs["TTG"]
    else:
        start_codons["TTG"][row] = 0

# Computing uniformity of start-codon per ortologus row
start_codons["uniformity"] = "NA"
for row in range(len(start_codons)):
    if len(set(start_codons.iloc[row, 1])) == 1:
        start_codons.iloc[row, 5] = "uniform"
    else:
        start_codons.iloc[row, 5] = "non-uniform"

# Constructing table withh all data about the ortologus rows (including pararows)
start_codons2 = start_codons.merge(orto_rows, on="ortologus_row", how="outer")

# Creating a table, combining data about gene and its start-codone
strain_gene_row = start_codons2[["Strains", "Genes", "ortologus_row", "uniformity", "ATG", "GTG", "TTG"]]   # Other codons deleted
summary_rows = df1.merge(strain_gene_row, on="ortologus_row")

# Identifying non-uniform rows
# We need to align only rows with different start-codons.

non_uniform_or_list = list(set(start_codons2.query("uniformity == 'non-uniform'").ortologus_row))
for orto_row in tqdm(non_uniform_or_list):
    subset = summary_rows.query("ortologus_row == @orto_row")
    with open (f"../{folder_name}/data/multialignments/{folder_name}_withoutp_{str(orto_row)}.fasta", "w") as nucleotide_fasta:
        for index, row in subset.iterrows():
            if row['type_of_the_gene'] != "pseudogene":
                nucleotide_fasta.write(">")
                nucleotide_fasta.write(row["id"])
                nucleotide_fasta.write("_")
                nucleotide_fasta.write(row["source"])
                nucleotide_fasta.write("_")
                nucleotide_fasta.write(row["start_codone"])
                nucleotide_fasta.write("\n")
                nucleotide_fasta.write(row["n_sequence"])
                nucleotide_fasta.write("\n")


# #### Saving summary-rows and start-codons2

summary_rows.to_csv(f"../{folder_name}/data/summary_rows_prokka.csv", index=False)
start_codons2.to_csv(f"../{folder_name}/data/start_codons2_prokka.csv", index=False)

