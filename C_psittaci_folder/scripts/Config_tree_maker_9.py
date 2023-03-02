#!/usr/bin/env python
# coding: utf-8

import sys
import os
import re
from tqdm import tqdm

current_path = "../data/multialignments/"
multialigns = [align for align in os.listdir(current_path) if align[-3:] == "afa"]

pattern = r">[\w]+"
for align in tqdm(multialigns):
    with open(f"../data/configs/{align[:-4]}.txt", "w") as writed_align:
        writed_align.write("DATASET_BINARY\n")
        writed_align.write("SEPARATOR COMMA\n")
        writed_align.write("DATASET_LABEL,start codons\n")
        writed_align.write("FIELD_SHAPES,2,2,2\n")
        writed_align.write("FIELD_LABELS,ATG,GTG,TTG\n")
        writed_align.write("FIELD_COLORS,#ff0000,#00ff00,#0000ff\n")
        writed_align.write("DATA\n")
        with open(os.path.join(current_path, align), "r") as readed_align:
                texts = readed_align.read()
                result = re.findall(pattern, texts)
                for header in result:
                    writed_align.write(header[1:])
                    if header[-3:] == "ATG":
                        writed_align.write(",1,-1,-1\n")
                    elif header[-3:] == "GTG":
                        writed_align.write(",-1,1,-1\n")
                    elif header[-3:] == "TTG":
                        writed_align.write(",-1,-1,1\n")

