#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 07:44:33 2025

@author: jrwill
"""

import pandas as pd


dataroot = "/Users/jrwill/prog/Pytheas_Folder/Pytheas_Qt_root/Pytheas_Qt/pytheas_data/"

light_file = "rna_mod_defs_light.xlsx"
heavy_file = "rna_mod_defs_hyx_heavy.xlsx"

nt_file = "nucleotide_table.xlsx"
nt_df = pd.read_excel(dataroot + nt_file, "Sheet1")
nt_dict = nt_df.set_index("AMBER_3_letter_code").to_dict(orient = "index")



light_sheet_dict = pd.read_excel(dataroot + light_file, None)

# dict key is 3-letter code, dict value is df for that nt

hyx = ["N9", "N7", "N1", "N3"]

for key, df in light_sheet_dict.items():
    orig_base = nt_dict[key]["Originating_base"]
    if orig_base == "ADE" or orig_base == "GUA":  # purines
        for idx, row in df.iterrows():
            if row.Atom_name in hyx:
                # print(key, row)
                print(key, idx, light_sheet_dict[key].loc[idx,"Atom_symbol"])
                light_sheet_dict[key].loc[idx,"Atom_symbol"] = "N15"
                print(key, idx, light_sheet_dict[key].loc[idx,"Atom_symbol"])
                print()
            
with pd.ExcelWriter(dataroot + heavy_file) as writer:
        for nt, df in light_sheet_dict.items():
            df.to_excel(writer, nt, index = False)

