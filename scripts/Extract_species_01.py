#!/usr/bin/env python

# Importing packages
import string
import json
import time
import requests
import argparse
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
arguments = parser.parse_args()
taxon = arguments.taxon
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
        if specie.startswith("unclassified"):
            gpb_species.remove(specie)
            counter += 1
        elif specie.startswith("endosymbiont"):
            gpb_species.remove(specie)
            counter += 1
        elif specie.startswith("Candidatus"):
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

# Executing 
org_keys = list(gbp_taxids.keys())
n_threads = 24

with ThreadPoolExecutor(n_threads) as pool:                                
    results = pool.map(counting_assembles, org_keys) 
    print(type(results), list(results)[:5])

# Saving results
json_gbp_s = json.dumps(big_ass)
with open(f"../data/gpb_species.json", "w") as jsngbps:
    jsngbps.write(json_gbp_s)
