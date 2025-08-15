#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 15:31:51 2023

@author: jrwill
"""
from datetime import datetime


import os


from pytheas_global_vars import pgv

from pytheas_IO import read_pytheas_file

# from pytheas_modules.worksheet_functions import format_worksheet_columns, write_worksheet_rows
# from pytheas_modules.scoring_functions import ppm_range, ppm_offset, consecutive_series, sumI, sumI_all, n_calc, L_calc
# from pytheas_modules.digest import digest
# import pytheas_modules.digest_functions as dg
# import match_functions as ma
# from discovery_functions_new import *
# from pytheas_modules.setup import *
# from permutations_dfs import score_distribution_plot

# from discovery_functions_new import build_mass_dict, discover_spectra
from discovery_functions import (build_mass_dict, discover_spectra,
                                    unpack_master_match_dict)
from match_functions import output_match_dict_file, consolidated_match_output

# import pytheas_discovery_dev.Combination_sum_class_work
# import pytheas_discovery_dev.Combination_Product_work

# def get_validated_target(ms2_key):
#     if ms2_key in pgv.valid_dict.keys():
#         target =  pgv.valid_dict[ms2_key]["frag"]
#         end3, tseq, end5 = target.split("_")
#         target3 = dg.parse_mod_seq(tseq)
#     else:
#         target = "I_unknown_am"
#         target3 = ["IAM", "UNK", "KNO", "WN "]
#     return target, target3


def discovery():
    global sc_ctr, length_ctr
    
    t1 = datetime.now()

    read_pytheas_file("nt_def_file")
    # read_pytheas_file("modification_set_file")
    
    # print("discovery: set_dict", pgv.set_dict)
    # ms2_dict, ms2_key_dict = ma.read_MS_data_file(pgv.MS_data_file)  # ms2_key_dict is suitable for .json
    # pgv.ms2_dict = ms2_dict

    read_pytheas_file("MS_data_file")  # read in MS data to pgv.ms2_dict
    
    # print("pgv: ", pgv.__dict__.keys())
    
    pgv.n_ms2_keys = len(list(pgv.ms2_dict.keys()))
    print("number of spectra = ", pgv.n_ms2_keys)

    # pgv.modification_set = ["training_set"]
    # pgv.widget_dict["modification_set"].load_value([pgv.modification_set])
    
#     # if pgv.Xcorr =='y':
#     #     ma.fft_MS2()

# #TODO    # need to loop thru labels
    mass_dict = {key:val.mass for key,val in pgv.nt_fragment_dict["light"].items() if "end" not in pgv.nt_def_dict[key]["Type"]}
    
#TODO make option panel for mod selection  -- presets and checkbox
    print("discovery: modification_set", pgv.modification_set, type(pgv.modification_set))
    # for s in pgv.modification_set:
    #     print()
    #     print("discovery set ", s, type(s))
    #     # print(pgv.set_dict)
    #     print(pgv.set_dict[s])
    print("discovery: set dict", [pgv.set_dict[s] for s in pgv.modification_set])
    # base_list = [b for set_name in pgv.modification_set for b in pgv.set_dict[set_name]]
    base_list = [b for s in  pgv.modification_set for b in pgv.set_dict[s] ]
    print("mod list", base_list)
    pgv.iso_mass_dict, pgv.iso_mass_list, pgv.mod_mass_list = build_mass_dict(mass_dict, base_list)

#     # pgv.iso_mass_dict, pgv.iso_mass_list, pgv.mod_mass_list = build_mass_dict(mass_dict, pgv.natural + pgv.training_mods)
#     # pgv.iso_mass_dict, pgv.iso_mass_list, pgv.mod_mass_list = build_mass_dict(mass_dict, pgv.natural + pgv.trna_mods)
    print("iso_mass_dict", pgv.iso_mass_dict)
    print("iso_mass_list", pgv.iso_mass_list)
    print("mod_mass_list", pgv.mod_mass_list)
    
#     #TODO see if this actually helps in discovery mode....
    
#     tt = datetime.now()
    pgv.master_precursor_dict = {}  # temporary dict to accumulate precursors and avoid recalculation ()
        
    pgv.master_match_dict = discover_spectra()
    
#     print("discovery took ", datetime.now() - tt )
    pgv.unpacked_match_dict = unpack_master_match_dict(pgv.master_match_dict)
    
    # for midx, mdict in pgv.unpacked_match_dict.items(): # reformat lists and dicts in
    #     ma.format_match_entry(mdict)

    output_match_dict_file("discovery_output")

#     ma.match_output_keys(pgv.unpacked_match_dict)

    pgv.top_match_dict, pgv.seq_match_dict, pgv.match_dict = consolidated_match_output("consolidated_discovery_output")
 
    
#     pgv.write_json(pgv.master_match_dict, os.path.join(pgv.pytheas_data_folder, "master_match_dict.json"))
#     pgv.write_json(pgv.unpacked_match_dict, os.path.join(pgv.pytheas_data_folder, "unpacked_match_dict.json"))
#     pgv.write_json(pgv.match_dict, os.path.join(pgv.pytheas_data_folder, "match.json"))
#     # pgv.write_json(top_match_dict, os.path.join(pgv.working_dir, "top_match.json"))
#     # pgv.write_json(seq_match_dict, os.path.join(pgv.working_dir, "seq_match.json"))
#     pgv.write_json(pgv.match_dict, os.path.join(pgv.pytheas_data_folder, "discovery_match.json"))
#     pgv.write_json(pgv.iso_mass_dict, os.path.join(pgv.pytheas_data_folder, "iso_mass_dict.json"))

#     if pgv.plot_sequence_map == 'y':
#         max_seq_len = max([len(mdict["seq_list"]) for mdict in pgv.mol_dict.values()])
#         if max_seq_len < 100:
#             ma.make_sequence_plot("discovery_sequence_map")
#         else:
#             ma.make_long_sequence_plot("discovery_sequence_map")

#     print()
#     print("discovery took :", datetime.now() - t1)

    return

