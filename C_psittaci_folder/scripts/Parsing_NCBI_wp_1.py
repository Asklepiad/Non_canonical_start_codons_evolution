#!/usr/bin/env python:
# coding: utf-8

# Installing packages
#get_ipython().system('pip3 install pandas')
#get_ipython().system('pip3 install biopython')

# Importing packages
import Bio
import math
import os
import sys
import argparse
import shutil
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from itertools import chain
from tqdm import tqdm
import json

# Asserting email for Entrez

parser = argparse.ArgumentParser()

parser.add_argument("organism_name", type=str)
parser.add_argument("mail", type=str)

arguments = parser.parse_args()
 
organism_name = arguments.organism_name
Entrez.email = arguments.mail



def extract_insdc(links): 
    linkset = [ls for ls in links[0]['LinkSetDb'] if
              ls['LinkName'] == 'assembly_nuccore_insdc']
    if 0 != len(linkset):
        uids = [link['Id'] for link in linkset[0]['Link']]
    else:
        uids = 0
    return uids

# Creating a variables for futher work with links
organism = organism_name
db_search = "assembly"
db_current = "nucleotide"
print("variables_ok")

# First searching helps to count complete number of links
search_handle = Entrez.esearch(db=db_search, term=organism)
search_record = Entrez.read(search_handle)
count = int(search_record["Count"])
print("search1_ok")

# Second searching
search_handle = Entrez.esearch(db=db_search, term=organism, retmax=count)
search_record = Entrez.read(search_handle)
print("search2_ok")

# Getting summary about links
idlist = search_record["IdList"]
@@ -73,6 +77,7 @@ def extract_insdc(links):
    record = Entrez.read(handle)
    if record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus'] == "Complete Genome":
        complete_ids.append(ids)
print("summaty ok")

# Taking ids for fetching. It collected all non-duplicated links in nucleotide database from assembly database.
links = []
links_checked = []
n = 0
for complete_id in tqdm(complete_ids):
    link_handle = Entrez.elink(dbfrom=db_search, db=db_current, from_uid=complete_id)
    link_record = Entrez.read(link_handle)
    uids = extract_insdc(link_record)
    if uids != 0:
        for uid in uids:
            if uid not in links_checked:    # Checking for duplicates
                links_checked.append(uid)
                links.append((uid, n))
                cumulative = 1
            else:
                cumulative = 0
        n += cumulative
print("linking_ok")

# Collecting data about assemblies
gb_records = []
for link in tqdm(links):
    gb_handle = Entrez.efetch(db=db_current, rettype="gb", retmode="text", id=link[0])
    gb_record = SeqIO.read(gb_handle, 'genbank')
    gb_records.append((gb_record, link[1]))
print("fetching_ok")

dna_type, tuples, source_list = [], [], [] # Creating list for identyfing the number of every assemblie DNA molecules (chromosome and any plasmids)
orglist = re.split(" |_", organism)
if len(orglist[-1]) < 10:
    last_letter = len(orglist[-1])
else:
    last_letter = 9
name = f"{orglist[0][0]}_{orglist[-1][:last_letter]}"
for rec in tqdm(gb_records):
    source = rec[1] # Number of assembly
    if "plasmid" in rec[0].description:
        dna_type.append("plasmid")
    else:
        dna_type.append("chromosome")
    source_list.append(source) 
    number = source_list.count(source) # Counting the DNA molecule of assembly
    tuples.append((source, number))
    mask = f"{name}{source}_{number}"
    with open (f"../{name}/data/for_prokka_fasta/{mask}.fasta", "w") as for_prokka_fasta:  # It is firstly needed to have a directory
        for_prokka_fasta.write(">")
        for_prokka_fasta.write(mask)
        for_prokka_fasta.write("\n")
        for_prokka_fasta.write(str(rec[0].seq))
        for_prokka_fasta.write("\n")

# Saving plasmid_code
jsonpc1 = json.dumps(tuples)
with open(f"../{name}/data/{name}_plasmid_code1.json", "w") as jspc1:
    jspc1.write(jsonpc1)
jsonpc2 = json.dumps(dna_type)
with open(f"../{name}/data/{name}_plasmid_code2.json", "w") as jspc2:
    jspc2.write(jsonpc2)
