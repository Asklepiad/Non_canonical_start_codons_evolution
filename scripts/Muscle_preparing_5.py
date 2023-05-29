#!/usr/bin/env python
# coding: utf-8

# On the next stage I need to download file with ortologus rows in the system. I have mounted google drive with project files (in the future I'm going to use GitHub for these purposes). 

import os
import sys
import re
import time
import argparse
import shutil
import json
import numpy as np
import pandas as pd
from tqdm import tqdm

# Preparing
start = time.time()

parser = argparse.ArgumentParser()

parser.add_argument("json_path", type=str)

arguments = parser.parse_args()

json_path = arguments.json_path
with open(os.path.join("../data/jsons", json_path), "r") as json_organism:
    json_organism = json.load(json_organism)

organism_name = json_organism[0]

lst = re.split(" |_", organism_name)
if len(lst[-1]) < 10:
    last_letter = len(lst[-1])
else:
    last_letter = 9
folder_name = f"{lst[0][0]}_{lst[-1][:last_letter]}"

preparing_finish = time.time()
print(f"Preparing time = {preparing_finish - start}")

# ### Orto rows with singletons, but without paralog-containing rows
# Uploading and modificating orto rows tables
orto_rows = pd.read_csv(f"../{folder_name}/data/{folder_name}.proteinortho.tsv", sep="\t") # Uploading dataframe with orto rows
orto_rows = orto_rows.rename(columns = {"# Species": "Species"})
orto_rows["ortologus_row"] = orto_rows.index + 1 # Creating ortorows numbers
orto_number = time.time()
print(f"Ortorow number creation = {orto_number - preparing_finish}")

# Uploading first csv-table and creating a new column in it
orto_rows_list = orto_rows.index
df1 = pd.read_csv(f"../{folder_name}/data/First_table.csv")

# Filling orto_row column (sounds like an oxymoron)
assemblies = orto_rows.columns.drop(['Species', 'Genes', 'Alg.-Conn.', 'ortologus_row'])
row_assigner = pd.melt(orto_rows, id_vars=["ortologus_row"], value_vars=assemblies)
row_assigner = row_assigner.query("value != '*'")
row_assigner["value"] = row_assigner["value"].str[:14]
row_assigner
row_assigner.drop(["variable"], axis="columns", inplace=True)
row_assigner.columns = ["ortologus_row", "id"]
row_assigner
df1 = df1.merge(row_assigner, on="id")
sec_orto_number = time.time()
print(f"Second ortorow number creation = {sec_orto_number - orto_number}")

# Number of the paralogs
print("Number of the paralogs =", sum(orto_rows.query("Genes > Species").Genes) - sum(orto_rows.query("Genes > Species").Species))

# Creating a subset without pararows
pararows_numbers = orto_rows.query("Genes > Species").ortologus_row
df1 = df1.query("ortologus_row not in @pararows_numbers").query("ortologus_row != 0")
paralogs_deleting = time.time()
print(f"Paralogs deleting = {paralogs_deleting - sec_orto_number}")

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
ortorow_sc = time.time()
print(f"Start-codons of ortologus row detecting = {ortorow_sc - paralogs_deleting}")



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
sc_freqs = time.time()
print(f"Computing of start-codons in each ortologus row = {sc_freqs - ortorow_sc}")

# Computing uniformity of start-codon per ortologus row
uniform = []
start_codons.start_codons.apply(lambda x: uniform.append("different") if len(x) > 1 else uniform.append("same"))
uniformity = pd.Series(uniform)
start_codons["uniformity"] = uniformity
print(f"Computing of uniformity = {unif_comp - sc_freqs}")

# Constructing table withh all data about the ortologus rows (including pararows)
start_codons2 = start_codons.merge(orto_rows, on="ortologus_row", how="outer")
first_merging = time.time()
print(f"First merging = {first_merging - sc_freqs}")

# Creating a table, combining data about gene and its start-codone
strain_gene_row = start_codons2[["Species", "Genes", "ortologus_row", "uniformity", "ATG", "GTG", "TTG"]]   # Other codons deleted
summary_rows = df1.merge(strain_gene_row, on="ortologus_row")
second_merging = time.time()
print(f"Second merging = {second_merging - first_merging}")

# Assigning cog to orto_row
row_cog = row_assigner.merge(df1[["id", "cog"]], on="id").groupby("ortologus_row").agg({"cog": "max"})
start_codons2 = start_codons2.merge(row_cog, on="ortologus_row")

# Adding organism_name to datasets
summary_rows["organism"] = folder_name
start_codons2["organism"] = folder_name

# Creating dataframe about ortorows with consistent number of columns
uniq_cog = start_codons2[["ortologus_row", "ATG", "GTG", "TTG", "uniformity", "Species", "Genes", "cog", "organism"]]

# We need to align only rows with different start-codons.

non_uniform_or_list = list(set(start_codons2.query("uniformity == 'different'").ortologus_row))
for orto_row in tqdm(non_uniform_or_list):
    subset = summary_rows.query("ortologus_row == @orto_row")
    with open (f"../{folder_name}/data/multialignments/{folder_name}_withoutp_{str(orto_row)}.fasta", "w") as nucleotide_fasta:
        for index, row in subset.iterrows():
            #if row['type_of_the_gene'] != "pseudogene":
            nucleotide_fasta.write(">")
            nucleotide_fasta.write(row["id"])
            nucleotide_fasta.write("_")
            nucleotide_fasta.write(row["source"])
            nucleotide_fasta.write("_")
            nucleotide_fasta.write(row["start_codone"])
            nucleotide_fasta.write("\n")
            nucleotide_fasta.write(row["n_sequence"])
            nucleotide_fasta.write("\n")
multial_creating = time.time()
print(f"Multialignment creating = {multial_creating - second_merging}")

# #### Saving summary-rows and start-codons2

summary_rows.to_csv(f"../{folder_name}/data/summary_rows_prokka.csv", index=False)
start_codons2.to_csv(f"../{folder_name}/data/start_codons2_prokka.csv", index=False)
uniq_cog.to_csv(f"../{folder_name}/data/{folder_name}_uniq_cog.csv", index=False)
saving_datasets = time.time()
print(f"Saving datsets = {saving_datasets - multial_creating}")
