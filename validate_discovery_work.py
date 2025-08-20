#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 14:46:29 2025

@author: jrwill
"""


m_df = upm_df[upm_df["ms2_key"] == 12]

print(max(test_df["Sp"]))

print(test_df["Sp"])

from pytheas_IO import read_json
import os
from digest_functions import generate_mod_seq

load_dir = "/Users/jrwill/prog/Pytheas_Folder/Pytheas_Data/Example_data/Discovery_Job_599"
jfile = "unpacked_match.json"
dis_dict = read_json(os.path.join(load_dir,"pytheas_json_files", jfile))

dis_df = pd.DataFrame.from_dict(dis_dict, orient = "index")


for idx, row in test_df.iterrows():
    # seq3 = row["frag3"][1:-2]
    # print(seq3)
    # seq3s = [str(s) for s in seq3]
    print(generate_mod_seq(row["frag3"]), row["ms2_key"], row["Sp"])
    
d_df = dis_df[dis_df["ms2_key"] == 12]

for idx, row in d_df.iterrows():
    # seq3 = row["frag3"][1:-2]
    # print(seq3)
    # seq3s = [str(s) for s in seq3]
    print(generate_mod_seq(row["frag3"]), row["ms2_key"], row["Sp"])

# %%

upm_df = pd.DataFrame.from_dict(pgv.unpacked_match_dict, orient = "index") # match output
dis_df = pd.DataFrame.from_dict(dis_dict, orient = "index") # discovery output



ms2_keys = list(set(list(upm_df["ms2_key"].unique()) + list(dis_df["ms2_key"].unique())))

n_matches = 0
n_sub_matches = 0

for key in ms2_keys:
    m_df = upm_df[upm_df["ms2_key"] == key]
    d_df = dis_df[dis_df["ms2_key"] == key]
    
    if len(d_df) != 0:
        d_top = d_df.iloc[0]
        d_seq = generate_mod_seq(d_top["frag3"])
     
    if len(m_df) != 0:
        m_top = m_df.iloc[0]
        m_seq = generate_mod_seq(m_top["frag3"])

    if len(d_df) == 0:
        print("ms2_key ", key, "no discovery", "match = ", m_seq)
        continue
    
    if len(m_df) == 0:
        print("ms_key ", key, "no match", "discovery = ", d_seq)
        continue
    
    
    if m_seq == d_seq:
        print("ms2_key", key, "match = discovery", m_seq)
        n_matches += 1
    else:
        print("ms2_key", key, "match != discovery ", m_seq, d_seq)
        ctr = 0
        for idx, row in d_df.iterrows():
            d_seq = generate_mod_seq(row["frag3"])
            ctr += 1
            if d_seq == m_seq:
                print("found match at row ", ctr)
                n_sub_matches += 1
                

print()
print("total matches = ", n_matches)
print("suboptimal matches ", n_sub_matches)
      
# %%

def validate_discovery():
    upm_df = pd.DataFrame.from_dict(pgv.unpacked_match_dict, orient = "index") # match output
    dis_df = pd.DataFrame.from_dict(pgv.unpacked_discovery_dict, orient = "index") # discovery output



    ms2_keys = list(set(list(upm_df["ms2_key"].unique()) + list(dis_df["ms2_key"].unique())))

    n_matches = 0
    n_sub_matches = 0

    for key in ms2_keys:
        m_df = upm_df[upm_df["ms2_key"] == key]
        d_df = dis_df[dis_df["ms2_key"] == key]
        
        if len(d_df) != 0:
            d_top = d_df.iloc[0]
            d_seq = generate_mod_seq(d_top["frag3"])
         
        if len(m_df) != 0:
            m_top = m_df.iloc[0]
            m_seq = generate_mod_seq(m_top["frag3"])

        if len(d_df) == 0:
            print("ms2_key ", key, "no discovery", "match = ", m_seq)
            continue
        
        if len(m_df) == 0:
            print("ms_key ", key, "no match", "discovery = ", d_seq)
            continue
        
        
        if m_seq == d_seq:
            print("ms2_key", key, "match = discovery", m_seq)
            n_matches += 1
        else:
            print("ms2_key", key, "match != discovery ", m_seq, d_seq)
            ctr = 0
            for idx, row in d_df.iterrows():
                d_seq = generate_mod_seq(row["frag3"])
                ctr += 1
                if d_seq == m_seq:
                    print("found match at row ", ctr)
                    n_sub_matches += 1
                    

    print()
    print("total matches = ", n_matches)
    print("suboptimal matches ", n_sub_matches)
          