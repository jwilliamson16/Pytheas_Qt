#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 15:31:51 2023

@author: jrwill
"""
from datetime import datetime
import os
from pathlib import Path

from pytheas_global_vars import pgv, pgc
from pytheas_IO import read_pytheas_file, save_json_files, load_pickle

from discovery_functions import (build_mass_dict, discover_spectra,
                                    unpack_master_discovery_dict, 
                                    validate_discovery)
from match_functions import output_match_dict_file, consolidated_match_output


def discovery():

    start = datetime.now()
    read_pytheas_file("MS_data_file")  # read in MS data to pgv.ms2_dict
    
    pgv.n_ms2_keys = len(list(pgv.ms2_dict.keys()))
    print("number of spectra = ", pgv.n_ms2_keys)

# TODO... implement Xcorr?    
#     # if pgv.Xcorr =='y':
#     #     ma.fft_MS2()
    
#TODO make option panel for mod selection  -- presets and checkbox
    # print("discovery: modification_set", pgv.modification_set, type(pgv.modification_set))
    # base_list = [b for s in  pgv.modification_set for b in pgv.set_dict[s] ]
    # print("mod list", base_list)
    
    # # #TODO    # need to loop thru labels
    # # TODO move to build_mass_dict??

#TODO fix 7MG and more generally charge on nucleotide
    build_mass_dict()
    
    pgv.ms2_file_list = discover_spectra()
    
    pgv.discovery_dict = {}
    for ms2_file in pgv.ms2_file_list:
        ms2_key = int(ms2_file.split("/")[-1].split(".")[0].split("_")[-1])
        pgv.discovery_dict[ms2_key] = load_pickle(ms2_file)

    pgv.unpacked_discovery_dict = unpack_master_discovery_dict(pgv.discovery_dict)
    
    discovery_job = pgv.discovery_job.split("_")[-1]
    discovery_file = "discovery_output" + "_" + discovery_job
    output_match_dict_file(pgv.unpacked_discovery_dict, discovery_file)

    pgv.top_discovery_dict, pgv.seq_discovery_dict, pgv.discovery_dict = consolidated_match_output(pgv.unpacked_discovery_dict, "consolidated_discovery_output")
 
#TODO  add Sp to validation output
    validate_discovery()
    
    json_dir = os.path.join(pgv.job_dir, "pytheas_json_files")
    Path(json_dir).mkdir(parents=True, exist_ok=True)
 
    save_json_files(pgc.discovery_json, json_dir)

# TODO try to match on sequence
#     if pgv.plot_sequence_map == 'y':
#         max_seq_len = max([len(mdict["seq_list"]) for mdict in pgv.mol_dict.values()])
#         if max_seq_len < 100:
#             ma.make_sequence_plot("discovery_sequence_map")
#         else:
#             ma.make_long_sequence_plot("discovery_sequence_map")

    print()
    print("discovery took :", datetime.now() - start)

    return

