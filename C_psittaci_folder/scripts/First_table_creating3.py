#!/usr/bin/env python
# coding: utf-8

# ### Uploading the data

import os
import sys
import argparse
import shutil
import json
import re
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio.Seq import Seq
from itertools import chain
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()

parser.add_argument("folder_name", type=str)
#parser.add_argument("cog_source", type=str, default="https://www.ncbi.nlm.nih.gov/research/cog/cogcategory/")

arguments = parser.parse_args()

folder_name = arguments.folder_name
#cog_url = arguments.cog_source

# Function for parsing assembly number
def get_assembly_number(pre_pattern1, string):
    pre_pattern2 = r"_[\d]+?\.gbk"
    first_bound = re.search(pre_pattern1, string)
    second_bound = re.search(pre_pattern2, string)
    result = string[first_bound.end():(second_bound.start())]
    return result


# Uploading the data 
current_path = f"../{folder_name}/data/{folder_name}_annotate/" # Path directory
pre_pattern1 = f"{folder_name}"
gb_records= []  # Reasserting the gb_records list
gbk_files = os.listdir(current_path)

for gbk in tqdm(gbk_files):
    assembly_number = get_assembly_number(pre_pattern1, gbk)
    link = f"{current_path}{gbk}"
    gbk_object = SeqIO.read(link, 'genbank')
    gb_records.append((gbk_object, assembly_number))

# Uploading info about which DNA is plasmid
with open(f"../{folder_name}/data//{folder_name}_plasmid_code1.json", "r") as tuples:
    tuples = json.load(tuples)
with open(f"../{folder_name}/data/{folder_name}_plasmid_code2.json", "r") as dna_type:
    dna_type = json.load(dna_type)
tuples = [tuple(x) for x in tuples]
plasmid_code = dict(zip(tuples, dna_type))

# Parsing
locus_tag, start_codone, source, n_sequence, aa_sequence, DNA_source, \
product, cog, p_c_unity, first_n, last_n, strand_n, \
length_n = [], [], [], [], [], [], [], [], [], [], [], [], []
for record_number in tqdm(range(len(gb_records))):
    for feature_number in range(1, len(gb_records[record_number][0].features)):
        if gb_records[record_number][0].features[feature_number].type == "CDS":
            
            ## Column 1 -- uniq id
            if 'locus_tag' in gb_records[record_number][0].features[feature_number].qualifiers.keys():
                locus_tag.append(gb_records[record_number][0].features[feature_number].qualifiers['locus_tag'][0])
            else:
                locus_tag.append("absent")
                
            ## Column 2 -- source
            source.append(gb_records[record_number][0].name)
            
            ## Column 3 -- nucleotide sequence & 7 -- start-codone
            first = gb_records[record_number][0].features[feature_number].location.start
            last = gb_records[record_number][0].features[feature_number].location.end
            strand = gb_records[record_number][0].features[feature_number].location.strand
            first_n.append(first)
            last_n.append(last)
            strand_n.append(strand)
            length_n.append(last - first)
            if strand == 1:
                sequence = str(gb_records[record_number][0].seq[first:last])
                n_sequence.append(sequence)
                start_codone.append(sequence[0:3])
            elif strand == -1:
                sequence = str(gb_records[record_number][0].seq.complement()[first:last])[::-1]
                n_sequence.append(sequence)
                start_codone.append(sequence[0:3])
            else:
                n_sequence.append("absent")
                
            ## Column 4 -- aminoacid sequence
            if 'translation' in gb_records[record_number][0].features[feature_number].qualifiers.keys():
                aa_sequence.append(gb_records[record_number][0].features[feature_number].qualifiers['translation'][0])
            else:
                aa_sequence.append("absent")
            
            ## Column 5 -- plasmid or chromosome
            index = re.search(r"[\d]+?_[\d]+$", gb_records[record_number][0].name)
            pc_key = tuple([int(i) for i in index.group().split("_")])
            if plasmid_code[pc_key] == "plasmid":
                DNA_source.append("plasmid")
            else:
                DNA_source.append("chromosome")
         
            # Column 6 -- protein, coding by this gene
            if 'product' in gb_records[record_number][0].features[feature_number].qualifiers.keys():
                product.append(gb_records[record_number][0].features[feature_number].qualifiers['product'][0])
            else:
                product.append("absent")
                
            ## Column 7 -- COG number
            if 'db_xref' in gb_records[record_number][0].features[feature_number].qualifiers.keys():
                cog.append(gb_records[record_number][0].features[feature_number].qualifiers["db_xref"][0][4:])
            else:
                cog.append("absent")
            
            ## Column 8 -- number of molecule (for identification and further plasmid marking)
            p_c_unity.append(gb_records[record_number][1])
            
df1 = pd.DataFrame(
    {
        "id": pd.Series(locus_tag),
        "source": pd.Series(source),
        "n_sequence": pd.Series(n_sequence),
        "aa_sequence": pd.Series(aa_sequence),
        "type_of_DNA_source": pd.Series(DNA_source),
        "start_codone": pd.Series(start_codone),
        "product": pd.Series(product),
        "cog": pd.Series(cog),
        "p_c_unity": pd.Series(p_c_unity),
        "first": pd.Series(first_n),
        "last": pd.Series(last_n),
        "strand": pd.Series(strand_n),
        "length": pd.Series(length_n),
    }
)


# For recoding COGs from ids to categories I have downloaded tsvs with all categories and appropriate COGs and have created python dictionary.
# Firstly I have uploaded tsvs to folder on google drive.

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


# On the next stage I have created 25 new columns for dataset, each contained data if protein has this COG.

# Adding 25 columns to dataset
for key in tqdm(cog_dict.keys()):   # Passing over COG categories
    cog_list = []             # For each category create list
    for cog_id in cog:        # Passing out every item
        if cog_id not in cog_dict[key]:     # For every gene, if gene's id not in dict
            cog_list.append(0)              # Mark it
        else:                               
            cog_list.append(1)              # Else mark another
    df1[key] = pd.Series(cog_list)   # Creating 25 new columns



df1.to_csv(f"../{folder_name}/data/First_table.csv", index=False)


for source in tqdm(set(df1["p_c_unity"])):
    subset = df1[df1["p_c_unity"] == source]
    with open (f"../{folder_name}/data/orto_rows/{folder_name}{str(source)}.fasta", "w") as protein_fasta:
        for index, row in subset.iterrows():
#            if row['type_of_the_gene'] != "pseudogene":
            protein_fasta.write(">")
            protein_fasta.write(row["id"])
            protein_fasta.write("_")
            protein_fasta.write(row["source"])
            protein_fasta.write("\n")
            protein_fasta.write(row["aa_sequence"])
            protein_fasta.write("\n")

