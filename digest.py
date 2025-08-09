#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 10:30:53 2023

@author: jrwill
"""


import os
from datetime import datetime
from pathlib import Path


from pytheas_global_vars import pgv, pgc
from digest_functions import (add_modifications, enzyme_digest, add_precursor_ions,
                            build_frag_dict, build_precursor_dict,  output_digest_file)
from pytheas_IO import read_pytheas_file, save_json_files

def digest():

    main_start = datetime.now()

    print()
    print("STEP 1: READ SEQUENCES")

    read_pytheas_file("fasta_file")
    print("   Read ", len(pgv.mol_dict.keys()), " raw sequences from ", pgv.fasta_file)
    print("   5' ends are ", pgv.mol_end5, "   3' ends are ", pgv.mol_end3)

    print()
    print("STEP 2: MODIFICATIONS")

    pgv.active_mods = add_modifications()
    if pgv.rna_mods == 'none':
        print("   No modifications present")
    else:
        print("   Modifications taken from ", pgv.rna_mods)
        print("  ", len(pgv.mod_dict.keys()), " modified positions found")
        print("  ", len(pgv.active_mods), " types of active modifications")
  
    print()
    print("STEP 3:  ENZYME DIGESTION")

    nfrags, pgv.unique_frag_dict = enzyme_digest()
    nprec = add_precursor_ions(pgv.unique_frag_dict, "all")

    if pgv.enzymes == ['none']:
        print("   No enzyme digestion performed")
    else:
        print("   Enzyme digestion performed with", pgv.enzymes)
    print("   ", len(pgv.unique_frag_dict.keys()), "unique fragments generated")
    print("   ",  nprec, " total precursor ions generated")

    print()
    print("STEP 4: NONREDUNDANT FRAGMENTS") 
      
    pgv.frag_dict = build_frag_dict(pgv.unique_frag_dict)
    n_frag_keys = len(pgv.frag_dict.keys())    
    
    print("   Number of nonredundant fragments in fragment dictionary: ", n_frag_keys)
    
    print()
    print("STEP 5:" "PRECURSOR ION GENERATION")
    
    pgv.unique_precursor_dict = build_precursor_dict()
 
    print("    Number of unique precursor ions: ", len(pgv.unique_precursor_dict))
    print()
    
    json_dir = os.path.join(pgv.job_dir, "pytheas_json_files")
    Path(json_dir).mkdir(parents=True, exist_ok=True)
 
    save_json_files(pgc.digest_json, json_dir)
    
    digest_job = pgv.digest_job.split("_")[-1]
    digest_file = "digest_output" + "_" + digest_job
    output_digest_file(digest_file)
    
    if pgv.run == "CL":
        pgv.widget_dict["digest_output"].load_value(digest_file + ".xlsx")
    else:
        pgv.digest_file = digest_file + ".xlsx"
    print("     In Silico Digest took: ", datetime.now() - main_start)
    return None
