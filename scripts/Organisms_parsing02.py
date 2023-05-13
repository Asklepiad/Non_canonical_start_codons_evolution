#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import shutil
import re
import json
from itertools import chain

parser = argparse.ArgumentParser()
parser.add_argument("json_path", type=str)
arguments = parser.parse_args()
json_path = os.path.abspath(arguments.json_path)


with open(json_path, "r") as organism_list:
    organism_list = json.load(organism_list)

for bac, ids in organism_list.items():
    orglist = re.split(" |_", bac)
    if len(orglist[-1]) < 10:
        last_letter = len(orglist[-1])
    else:
        last_letter = 9
    name = f"{orglist[0][0]}_{orglist[-1][:last_letter]}"
    if not os.path.exists(f"../data/jsons/{name}.json"):
        jdir = [bac, {bac: ids}]
        new_file = json.dumps(jdir)
        with open(f"../data/jsons/{name}.json", "w") as jspc2:
            jspc2.write(new_file)

