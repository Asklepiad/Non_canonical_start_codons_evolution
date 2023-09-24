#!/usr/bin/env python

# Importing packages
import string
import json
import time
import requests
import argparse
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from threading import Thread
from tqdm import tqdm
from ete3 import NCBITaxa
from Bio import Entrez
from bs4 import BeautifulSoup

# Preparing databases
Entrez.email = "bogdan.sotnikov.1999@mail.ru"
ncbi = NCBITaxa()

# Parsing argument
parser = argparse.ArgumentParser()
parser.add_argument("taxon", type=str)
parser.add_argument("sp_per_genus", type=int)
arguments = parser.parse_args()
taxon = arguments.taxon
sp_per_genus = arguments.sp_per_genus
print(f"Counting number of species in {taxon}...")

# Counting number of taxons
# Limit = species
descendants = ncbi.get_descendant_taxa(taxon, intermediate_nodes=True, collapse_subspecies=True)
tmp = ncbi.translate_to_names(descendants)
print(f"Number of taxons (the lowerest level - specie) of {taxon} is {len(tmp)}")

# Limit = genus
descendants2 = ncbi.get_descendant_taxa(taxon, intermediate_nodes=True, rank_limit="genus")
tmp2 = ncbi.translate_to_names(descendants2)
print(f"Number of taxons (the lowerest level - genus) of {taxon} is {len(tmp2)}")

# Creating list with all species, but not upper taxons
gpb_species = [specie for specie in tmp if specie not in tmp2]
print(f"Number of {taxon} species before sorting is {len(gpb_species)}")

# Delete uncultured, unclassified and endosymbiont bacteria
counter = 1
while counter > 0:
    counter = 0
    for specie in gpb_species:
        if specie.__contains__("unclassified"):
            gpb_species.remove(specie)
            counter += 1
        elif specie.__contains__("endosymbiont"):
            gpb_species.remove(specie)
            counter += 1
        elif specie.__contains__("Candidatus"):
            gpb_species.remove(specie)
            counter += 1
        elif specie.__contains__(" sp."):
            gpb_species.remove(specie)
            counter += 1
        elif specie.__contains__("uncultured"):
            gpb_species.remove(specie)
            counter += 1
print(f"Number of {taxon} species after first sorting is {len(gpb_species)}")


# Deleting species groups and subgroups
gbp_taxids = ncbi.get_name_translator(gpb_species)
deleting_taxa = []
for value in gbp_taxids.values():
    new_val = ncbi.get_rank(value)
    if 'species' not in new_val.values():
        for i in new_val.keys():
            deleting_taxa.append(i)

deleting_keys = []
for key, value in gbp_taxids.items():
    if value[0] in deleting_taxa:
        deleting_keys.append(key)
        
for key in deleting_keys:
    del gbp_taxids[key]
print(f"Number of {taxon} species after second sorting is {len(gbp_taxids)}")

# Preparing for multiprocessing execution
url = "https://www.ncbi.nlm.nih.gov/assembly"
del_ass = []
big_ass = {}
def counting_assembles(org_name):
    try:
        response = requests.get(url, params={'term': f'((("chromosome"[Assembly Level]) OR "complete genome"[Assembly Level]) OR "scaffold"[Assembly Level]) AND {org_name}[Organism]'})
        ass_soup = BeautifulSoup(response.content, "lxml")
        pre_complete_a = ass_soup.findChildren("a", attrs={"data-value_id": "complete"})[0]
        complete_a = int(pre_complete_a.find_next_sibling().text[1:-1].replace(",", ""))
        pre_chrom_a = ass_soup.findChildren("a", attrs={"data-value_id": "chromosome"})[0]
        chrom_a = int(pre_chrom_a.find_next_sibling().text[1:-1].replace(",", ""))
        pre_scaff_a = ass_soup.findChildren("a", attrs={"data-value_id": "scaffold"})[0]
        scaff_a = int(pre_scaff_a.find_next_sibling().text[1:-1].replace(",", ""))
        ass = complete_a + chrom_a + scaff_a
        print(f"Organism: {org_name}\nNumber of assemblies: {ass}")
        if ass < 4:
            del_ass.append(org_name)
        else:
            big_ass[org_name] = ass
    except (Exception) as err:
        print(f"Problem with {org_name}")
        print(err)

# Executing 
org_keys = list(gbp_taxids.keys())
n_threads = 24

with ThreadPoolExecutor(n_threads) as pool:                                
    results = pool.map(counting_assembles, org_keys) 
    print(type(results), list(results)[:5])

# Saving results
json_gbp_s = json.dumps(big_ass)
with open(f"../data/pre_gpb_species.json", "w") as jsngbps:
    jsngbps.write(json_gbp_s)

genuses = np.unique([specie.split()[0] for specie in big_ass.keys()])
big_result = {}
for genus in genuses:
    per_gen_dict = {}
    for specie in big_ass:
        if specie.startswith(genus):
            per_gen_dict[specie] = big_ass[specie]
    if len(per_gen_dict) >= sp_per_genus:
        for repeat in range(sp_per_genus):
            spec = max(per_gen_dict, key=per_gen_dict.get)
            big_result[spec] = per_gen_dict[spec]
            del per_gen_dict[spec]
    elif len(per_gen_dict) > 0:
        for repeat in range(len(per_gen_dict)):
            spec = max(per_gen_dict, key=per_gen_dict.get)
            big_result[spec] = per_gen_dict[spec]
print(f"Number of species for further analysis: {len(big_result)}")

json_gbp_s_new = json.dumps(big_result)
with open(f"../data/gpb_species.json", "w") as jsngbps_n:
    jsngbps_n.write(json_gbp_s_new)
