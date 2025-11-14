#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 20:20:41 2023

@author: jrwill

"""

# This is the build for the rna_mod_defs_light/heavy from original sources
import pandas as pd
import numpy as np
from pandas import ExcelWriter


dataroot = "/Users/jrwill/prog/pytheas_tk_interface/pytheasroot/pytheasdata/mod_nucleotide_data_prep/"
topfile = "all_modrna08.lib"

with open(dataroot + topfile) as file:
    lines = [line.rstrip() for line in file]
    
modnt_list = []
idx = 1
line =  lines[idx]

while "!" not in line: # index of nucleotides
    modnt_list.append(line.replace("\"","").strip())
    idx +=1
    line = lines[idx]

nt_dict = {}
topo_dict = {}
idx = 0
for line in lines: # find position of each nucleotide definition 
    if ".unit.atoms " in line:
        fields = line.split(".")
        nt = fields[1]
        nt_dict[nt] = idx
    idx += 1

# find complete set of atom names

atom_list = []
for nt in modnt_list:
    idx = nt_dict[nt] + 1
    line = lines[idx]
    topo_dict[nt] = []
    while "!" not in line:
        fields = line.split(" ")
        atom = fields[1].replace("\"", "")
        atom_list.append(atom)
        print(fields[1],len(fields[1]),atom)
        topo_dict[nt].append(atom)
        idx += 1
        line = lines[idx]

atom_set = set(atom_list)

group_dict = {}
group_dict["P"] = ["P", "O1P", "O2P"]
group_dict["O5'"] = ["O5'"]
group_dict["O3'"] = ["O3'"]
group_dict["R"] = ["C1'","C2'", "C3'","C4'", "O4'", "C5'", "H1'", "H2'", "O2'", "HO'2", "H3'", "H4'", "H5'1", "H5'2", "CM2", "HM'1", "HM'2", "HM'3"]
# if not these, it's "B"

mod_nt_file = dataroot + "amber_rna_mods_nts_JSL.xlsx"
mod_nt_df = pd.read_excel(mod_nt_file)

def atom_symbol(atom):
    asym = "X99"
    alett = [c for c in atom if c.isalpha()]
    if "C" in alett:
        asym = "C12"
    if "N" in alett:
        asym = "N14"
    if "P" in alett:
        asym = "P31"
    if "O" in alett:  # O trumps P
        asym = "O16"
    if "H" in alett:
        asym = "H1"   # last check for case of HO2'
    # print(atom, asym)
    return asym

topo_dict = {}
n_atoms = 1
for i, row in enumerate(mod_nt_df.itertuples(), 1):
    if pd.isna(row.Pytheas_1_letter_code):
        continue
    code3 = row.AMBER_3_letter_code
    code1 = row.Pytheas_1_letter_code
    print(i, row.Name, row.Pytheas_1_letter_code, code3, nt_dict[code3])
    idx = nt_dict[code3] + 1
    line = lines[idx]
    topo_dict[code3] = []
    while "!" not in line:
        fields = line.split(" ")
        atom = fields[1].replace("\"", "")
        # atom_list.append(atom)
        atom_group = "B"
        atom_symb = atom_symbol(atom)
        for key, alist in group_dict.items():
            if atom in alist:
                atom_group = key
        # print(fields[1],len(fields[1]),atom, atom_symb, atom_group)
        topo_dict[code3].append([atom, atom_symb, n_atoms, atom_group])
        if atom == "O2P" and atom_group == "P":
            topo_dict[code3].append(["HO2P","H1", 1, atom_group])  # add proton for neutral mass

        idx += 1
        line = lines[idx]

# convert each list into a dataframe
# df_dict = {nt: pd.DataFrame(data , columns=['Atom_name', 'Atom_symbol', 'Atom_group']) for nt, data in topo_dict.items() }    

# with ExcelWriter(dataroot + "rna_mod_defs.xlsx") as writer:
#         for nt, df in df_dict.items():
#             df.to_excel(writer, nt, index = False)
            
### repeat for four natural bases 

topfile = "all_nucleic94.lib" # AMBER file for natural nts

with open(dataroot + topfile) as file:
    lines = [line.rstrip() for line in file]
    
modnt_list = []
idx = 1
line =  lines[idx]

while "!" not in line: # index of nucleotides
    modnt_list.append(line.replace("\"","").strip())
    idx +=1
    line = lines[idx]

nt_dict = {}
idx = 0
for line in lines: # find position of each nucleotide definition 
    if ".unit.atoms " in line:
        fields = line.split(".")
        nt = fields[1]
        nt_dict[nt] = idx
    idx += 1

          
unmod_topo_dict = {}
for code3 in ["RA", "RC", "RG", "RU"]:
    idx = nt_dict[code3] + 1
    line = lines[idx]
    unmod_topo_dict[code3] = []
    while "!" not in line:
        fields = line.split(" ")
        atom = fields[1].replace("\"", "")
        # atom_list.append(atom)
        atom_group = "B"
        atom_symb = atom_symbol(atom)
        for key, alist in group_dict.items():
            if atom in alist:
                atom_group = key
        # print(fields[1],len(fields[1]),atom, atom_symb, atom_group)
        unmod_topo_dict[code3].append([atom,atom_symb, n_atoms, atom_group])
        if atom == "O2P" and atom_group == "P":
            unmod_topo_dict[code3].append(["HO2P","H1", 1, atom_group])  # add proton for neutral mass

        idx += 1
        line = lines[idx]

# PART 1:  add unmodified dataframes
df_dict = {}
name_dict = {"RA": "ADE", "RC": "CYT", "RG": "GUA", "RU": "URI"}

for nt in ["RA", "RC", "RG", "RU"]:
    newnt = name_dict[nt]
    df_dict[newnt] =  pd.DataFrame(unmod_topo_dict[nt] , columns=['Atom_name', 'Atom_symbol', 'N_atoms', 'Atom_group'])

# convert each list into a dataframe
# df_dict = {nt: pd.DataFrame(data , columns=['Atom_name', 'Atom_symbol', 'Atom_group']) for nt, data in topo_dict.items() }    

# PART 2: add modified dataframes

for nt, data in topo_dict.items():
    df_dict[nt] = pd.DataFrame(data , columns=['Atom_name', 'Atom_symbol', 'N_atoms', 'Atom_group'])

# PART 3: add patch residues

# read in patch dict

patch_file = dataroot + "patch_groups.xlsx"
patch_sheet_dict = pd.read_excel(patch_file, None)
for nt,nt_df in patch_sheet_dict.items():
    df_dict[nt] = nt_df

    
# PART 4: add isobaric tag residues
tag_file = dataroot + "isobaric_tag_groups.xlsx"
tag_sheet_dict = pd.read_excel(tag_file, None)
for nt,nt_df in tag_sheet_dict.items():
    df_dict[nt] = nt_df

  
with ExcelWriter(dataroot + "rna_mod_defs_light.xlsx") as writer:
        for nt, df in df_dict.items():
            df.to_excel(writer, nt, index = False)

# make heavy template with 15N
skip_list = ["TG1", "TG2", "TG3", "TG4"] # skip heavy 15N for isobaric tags

for key,df in df_dict.items():
    if key in skip_list:
        continue
    for idx,row in df.T.items():
        if row["Atom_symbol"] == "N14":
            row["Atom_symbol"] = "N15"
            # df_dict[key]["Atom_symbol"][idx] = "N15"
            # row["Atom_symbol"] = "N15"

with ExcelWriter(dataroot + "rna_mod_defs_heavy.xlsx") as writer:
        for nt, df in df_dict.items():
            df.to_excel(writer, nt, index = False)

# Have Excel file with all nucleotide topologies!!
                
    