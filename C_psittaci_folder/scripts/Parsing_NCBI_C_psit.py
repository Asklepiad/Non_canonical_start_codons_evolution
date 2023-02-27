#!/usr/bin/env python
# coding: utf-8

# # Bogdan Sotnikov

# Installing packages
get_ipython().system('pip3 install pandas')
get_ipython().system('pip3 install biopython')

# Importing packages
import Bio
import math
import os
import sys
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
from google.colab import files   # For working with colab only
from google.colab import drive   # For working with colab only

# Asserting email for Entrez
Entrez.email = "bogdan.sotnikov.1999@mail.ru"

# Google drive mounting
drive.mount('/content/gdrive')    # For working with colab only


# After reading documentation of [biopython tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf), [Entrez documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_) and some additional resources ([ref1](https://notebook.community/widdowquinn/Teaching-EMBL-Plant-Path-Genomics/worksheets/01-downloading_data_biopython), [ref2](https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python)) I had idea about the parsing pipeline.
# 
# 1. I need to use ESearch twice. Firstly for identifying the count of assembles of interest, and then with retmax equal count from the first search.
# The example below will be executed on the full-genome assembles of *Chlamydia psittaci* and after then will be rewrited to general form.

# I have known there are only 27 complete genome assemblies of *Chlamydia psittaci* and only 83 assemblies overall. But I have had 105 ids.
# 
# 2. The second step of pipeline was to identify, which of these ids belong to complete assemblied genomes. I used esummary for this purpose.
# 
# 3. On this stage I used [list of Entrez links' abbreviations](https://eutils.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html) and modified code [from this source](https://notebook.community/widdowquinn/Teaching-EMBL-Plant-Path-Genomics/worksheets/01-downloading_data_biopython) 
# 
# **NB: All th code between this block and day 4 is only of historical interest.**

# Final pipeline:
# 
# 1. Searching of *Chlamydia psittaci* ids.
# 2. Filtering only complete assemblies.
# 3. Creating a function for filtering only GenBank ids.
# 4. Taking accession numbers for GenBank and using function, created on the previous step.
# 5. Downloading .gbk files.
# 6. Parsing them and creating a pandas dataframe with locus_tag source of the gene, nucleotide and aminoacid sequence.

# Code from source above, with modifications. It is for choosing GenBank nucleotide data only.

def extract_insdc(links):
    '''
    Selects only links, refered to nuccore_insdc
        :param links: list of NCBI links, received by Elinks 
    :return uids: list of UIDs, located in nuccore
    '''
    linkset = [ls for ls in links[0]['LinkSetDb'] if
              ls['LinkName'] == 'assembly_nuccore_insdc']
    if 0 != len(linkset):
        uids = [link['Id'] for link in linkset[0]['Link']]
    else:
        uids = 0
    return uids

# Creating a variables for futher work with links
organism = "Chlamydia psittaci"
db_search = "assembly"
db_current = "nucleotide"

# First searching helps to count complete number of links
search_handle = Entrez.esearch(db=db_search, term=organism)
search_record = Entrez.read(search_handle)
count = int(search_record["Count"])

# Second searching
search_handle = Entrez.esearch(db=db_search, term=organism, retmax=count)
search_record = Entrez.read(search_handle)

# Getting summary about links
idlist = search_record["IdList"]
complete_ids = []   # List of full completed genomes' ids
for ids in tqdm(idlist):
    handle = Entrez.esummary(db=db_search, id=ids)
    record = Entrez.read(handle)
    if record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus'] == "Complete Genome":
        complete_ids.append(ids)

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

# Collecting data about assemblies
gb_records = []
for link in tqdm(links):
    gb_handle = Entrez.efetch(db=db_current, rettype="gb", retmode="text", id=link[0])
    gb_record = SeqIO.read(gb_handle, 'genbank')
    gb_records.append((gb_record, link[1]))

# Creating a dictionary, which makes matching berween number of record and type of DNA source: chromosome or plasmid
plasmid_code = {}
for record_number in range(len(gb_records)):
    if "plasmid" in gb_records[record_number][0].description:
        plasmid_code[record_number] = "plasmid"
    else:
        plasmid_code[record_number] = "chromosome"
