#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:42:23 2023

@author: jrwill
"""

from pypdf import PdfReader
import csv

reader = PdfReader(dataroot + "amber_rna_mod_nts_JSL.pdf")
text = ""
for page in reader.pages:
    text += page.extract_text() + "\n"

lines = text.split("\n")

idx = 0
for line in lines:
    if "Table 3" in line:
        print(idx)
    idx +=1
    
# Table 3 starts on line 144

# table body starts on line 148    

idx = 148
for line in lines[148:]:
    if "PDB reference is given" in line:
        print(idx)
    idx +=1

# table ends at line 214
dataroot = "/Users/jrwill/prog/pytheas_tk_interface/pytheasroot/pytheasdata/"

# write out temp file

tlines = lines[148:214]
with open(dataroot + "amber_rna_mods_nts_JSL.txt", "w") as csv_file:
    for line in tlines:      
        csv_file.write(line + "\n")
        
        
# had to manually edit file because many lines were concatenate
# .txt file has 107 entries = Table 3
# read back in for further patching

# read back in and capture 3-letter code
with open(dataroot + "amber_rna_mods_nts_JSL_fixed.txt") as file:
    tlines = [line.rstrip() for line in file]

tslines = []
for line in tlines:
    sline = line.split()
    lenlist = [len(s) for s in sline]
    print(sline,lenlist)
    if 3 in lenlist:
        threepos = lenlist.index(3)
        name = " ".join(sline[0:threepos])
        three = sline[threepos]
        
        tslines.append([name,three])  
    

with open(dataroot + "amber_rna_mods_nts_JSL.csv", "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for line in tslines:
            writer.writerow(line)

# This worked!


