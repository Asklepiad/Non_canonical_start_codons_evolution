#!/usr/bin/env python

# Importing packages
import os
import sys
import argparse
import re
import Bio
from Bio import SeqIO
from Bio import Entrez
from tqdm import tqdm
import json

# Asserting email for Entrez

parser = argparse.ArgumentParser()

parser.add_argument("path", type=str)
parser.add_argument("mail", type=str)

arguments = parser.parse_args()

json_file_name = arguments.path
Entrez.email = arguments.mail
json_path = f"../data/jsons/{json_file_name}.json"
with open(json_path, "r") as json_organism:
    json_organism = json.load(json_organism)
 
organism_name = json_organism[0]
pre_complete_ids = list(json_organism[1].values())[0]
if type(pre_complete_ids[len(pre_complete_ids)-1]) == list:
    scaffs = pre_complete_ids[len(pre_complete_ids)-1]
    del pre_complete_ids[len(pre_complete_ids)-1]
    pre_complete_ids += scaffs


def extract_insdc(links, db_search, complete_id):
    search_handle = Entrez.esummary(db=db_search, id=complete_id)
    search_record = Entrez.read(search_handle)
    ass_level = search_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']
    if ass_level in ["Chromosome", "Complete Genome"]:
        linkset = [ls for ls in links[0]['LinkSetDb'] if
              ls['LinkName'] == 'assembly_nuccore_insdc']
        if 0 != len(linkset):
            uids = [link['Id'] for link in linkset[0]['Link']]
        else:
            uids = 0
    elif ass_level == "Scaffold":
        linkset = [ls for ls in links[0]['LinkSetDb'] if
              ls['LinkName'] == 'assembly_nuccore_refseq']
        if 0 != len(linkset):
            uids = [link['Id'] for link in linkset[0]['Link']]
        else:
            uids = 0
    return uids

def download_links(db_search, db_current, complete_id, timer, num_link):
    if timer > 0:
        link_handle = Entrez.elink(dbfrom=db_search, db=db_current, from_uid=complete_id)
        link_record = Entrez.read(link_handle)
        uids = extract_insdc(link_record, db_search, complete_id)
        if uids != 0:
            for uid in uids:
                if uid not in links_checked:  # Checking for duplicates
                    links_checked.append(uid)
                    links.append((uid, num_link))
                    cumulative = 1
                else:
                    cumulative = 0
            num_link += cumulative
    return num_link

def forced_download_links(db_search, db_current, complete_id, timer, level, num_link):   # Наверное, это можно реализовать декоратором
    try:
        num_link = download_links(db_search=db_search, db_current=db_current, complete_id=complete_id, timer=timer, num_link=num_link)
    except RuntimeError:
        print(f"Problem is with {complete_id}")
        timer -= 1
        level += 1
        print(f"We are on the {level} level now")
        num_link = forced_download_links(db_search, db_current, complete_id, timer, level, num_link)
    return num_link

# Creating a variables for futher work with links
organism = organism_name
db_search = "assembly"
db_current = "nucleotide"
print("variables_ok")

if pre_complete_ids[0][0:3] == "GCF":
    complete_ids = []
    for gcf in pre_complete_ids:
        search_handle = Entrez.esearch(db_search, gcf)
        search_record = Entrez.read(search_handle)
        complete_ids.append(search_record["IdList"])
else:
    complete_ids = pre_complete_ids

# Taking ids for fetching. It collected all non-duplicated links in nucleotide databiase from assembly database.
# We use try-except for excepting problems with network temorary lags, which ruined our code
links = []
links_checked = []
num_link = 0
for complete_id in tqdm(complete_ids):
    timer = 5
    level = 1
    num_link = forced_download_links(db_search, db_current, complete_id, timer, level, num_link)
print("linking_ok")

# Collecting data about assemblies
gb_records = []
for link in tqdm(links):
    gb_handle = Entrez.efetch(db=db_current, rettype="gb", retmode="text", id=link[0])
    gb_record = SeqIO.read(gb_handle, 'genbank')
    gb_records.append((gb_record, link[1]))
print("fetching_ok")

dna_type, tuples, source_list = [], [], []  # Creating list for identyfing the number of every assemblie DNA molecules (chromosome and any plasmids)
orglist = re.split(" |_", organism)
if len(orglist[-1]) < 10:   # 10 means only first nine letters (or lesser) will be in shortened name of bacteria. The cause of this approach is Prokka's bug.
    last_letter = len(orglist[-1])
else:
    last_letter = 9
name = f"{orglist[0][0]}_{orglist[-1][:last_letter]}"

for rec in tqdm(gb_records):
    try:
        source = rec[1]  # Number of assembly
        if "plasmid" in rec[0].description:
            dna_type.append("plasmid")
        else:
            dna_type.append("chromosome")
        source_list.append(source) 
        number = source_list.count(source)  # Counting the DNA molecule of assembly
        tuples.append((source, number))
        mask = f"{name}{source}_{number}"
        with open(f"../{name}/data/for_prokka_fasta/{mask}.fasta", "w") as for_prokka_fasta:  # It is firstly needed to have a directory
            for_prokka_fasta.write(">")
            for_prokka_fasta.write(mask)
            for_prokka_fasta.write("\n")
            for_prokka_fasta.write(str(rec[0].seq))
            for_prokka_fasta.write("\n")
    except Bio.Seq.UndefinedSequenceError:
        print(f"Mistake is in {rec[0].name}")
        continue

# Saving plasmid_code
jsonpc1 = json.dumps(tuples)
with open(f"../{name}/data/{name}_plasmid_code1.json", "w") as jspc1:
    jspc1.write(jsonpc1)
jsonpc2 = json.dumps(dna_type)
with open(f"../{name}/data/{name}_plasmid_code2.json", "w") as jspc2:
    jspc2.write(jsonpc2)
