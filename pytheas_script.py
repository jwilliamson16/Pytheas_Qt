#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 12:02:10 2025

@author: jrwill
"""

import sys
import multiprocessing
from discovery_functions import discover_spectra, build_mass_dict
from pytheas_IO import read_pytheas_file
import os
from pytheas_global_vars import pgv, pgc

def discovery_fork(ms2_key):
    global pgv
    # from pytheas_global_vars import pgv, pgc
    # pgv = PGV
    spec = pgv.ms2_dict[ms2_key]
    # if pgv.max_spectra != "all":
    #    if ms2_ctr > int(pgv.max_spectra):
    #        break
    
    # pgv.spec_dir = os.path.join(pgv.job_dir, "discovery_data_" + str(ms2_key))
    # pgv.plot_dir = pgv.spec_dir
    # Path(pgv.spec_dir).mkdir(parents=True, exist_ok=True)

    # if ms2_ctr > 300:
    #     break
    ms2_ctr += 1
    t1 = datetime.now()
    cidx = 0 # composition index
    master_composition_dict = {}
    composition_dict = {}
    if pgv.ion_mode == "-":  
        zsign = -1
    else:
        zsign = 1
         
    z = abs(spec.z)

    mobs = spec.mz1*z - zsign * z * pgv.hmass  # m0 for mz1
    print()
    print("*******")
    print("starting: ", ms2_key, spec.mz1, spec.z, zsign, z, mobs)
    for end5 in pgv.frag_end5:
        for end3 in pgv.frag_end3:
            end5n = pgv.end_dict["end5"][end5]
            end3n = pgv.end_dict["end3"][end3]
            for label in pgv.isotopic_species:

                mobs_adj = mobs - (pgv.nt_fragment_dict[label][end5n].mass + pgv.nt_fragment_dict[label][end3n].mass)
    
                seq_comps, ppm_tol, tol = find_compositions(mobs, mobs_adj) # need to add label 
                print("Found ", len(seq_comps), " Compositions:")
                for comp in seq_comps:
                    print("   ", comp)
                # print("m0 tolerance: ", ms2_key, ppm_tol, round(tol, 5), "# comps = ", len(seq_comps))
                if len(seq_comps) > 25:
                    print("TOO MANY COMPOSITIONS...Skip")
                    continue                    # print(seq_comps)
                composition_dict[ms2_key] = {"n_comps": len(seq_comps), "seq_comps":seq_comps, "ppm_tol": ppm_tol, "tol":tol} # needed??? output separately?
                composition_dict[ms2_key]["prec_dict"] = {"mz1": spec.mz1, "z": spec.z, "m0": mobs, "end5": end5, "end3": end3, "end5_3": end5n, "end3_3": end3n, 
                              "label": label}
 
                for seq_comp in seq_comps:  # match_composition() function
                    nperms = number_permutations(seq_comp)
                    # if nperms > pgv.max_permutations:
                    if len(seq_comp) > 5:
                        pgv.check_dfs_matching = 'y'
                        print()
                        print("using check_dfs on ms2_key ", ms2_key, " of length ", str(len(seq_comp)), "with", nperms, "permutations")
                        print()
                        # continue
                    else:
                        pgv.check_dfs_matching = 'y'
                    f3 = fragment_sequence([end5n] + seq_comp + [end3n])
                    precursor_dict, top_sort, top_keys, top_frags, top_match = match_permutations_dfs(f3, ms2_key, label, cidx, len(seq_comps))
                    unpacked_precursor_dict = prune_precursor_dict(precursor_dict)
                    # insert here:  polish permutations
                    
                    # maybe do Sp histogram first
                    master_composition_dict[cidx] = unpacked_precursor_dict
                    cidx += 1
    
    #TODO only save if there are matches
    if len(master_composition_dict) > 0:
        pgv.spec_dir = os.path.join(pgv.job_dir, "discovery_data_" + str(ms2_key))
        pgv.plot_dir = pgv.spec_dir
        Path(pgv.spec_dir).mkdir(parents=True, exist_ok=True)

        pickle_file = os.path.join(pgv.spec_dir, "discovery_" + str(ms2_key) + ".pkl")
        ms2_file_list.append(pickle_file)
        print("Saved: ", pickle_file)

        save_pickle(master_composition_dict, pickle_file)
        
        # print("master_composition_dict", master_composition_dict.keys())
        # for key, mdict in master_composition_dict.items():
        #     print(key,mdict.keys() )
        # print(master_composition_dict)
        discovery_file = "discovery_" + str(ms2_key)
        # unpacked_discovery_dict = unpack_master_discovery_dict(master_composition_dict)
        unpacked_discovery_dict = unpack_single_discovery_dict(ms2_key, master_composition_dict)
        # print("unpacked_discovery_dict", unpacked_discovery_dict)
        
        # TODO fix this hack 
        temp = pgv.job_dir
        pgv.job_dir = pgv.spec_dir
        output_match_dict_file(unpacked_discovery_dict, discovery_file)
        pgv.job_dir = temp

        print("spectrum ", ms2_key, " took ", datetime.now() - t1 )
        print("*******")
        print()
        return pickle_file


def init_worker(PGV):
    global pgv
    pgv = PGV

from pytheas_Qt import *

print("script sys argv: ", sys.argv)
print("**********loading parameters from script")
Load_Global_Vars()
print("************running digest from script")
# inSilicoDigest()
# matchSpectra()

pgv.max_mods = 1
# pgv.ms2_key_list = list(pgv.ms2_dict.keys())
pgv.ms2_key_list = list(range(10))

pgv.match_slope = 1.0 # need to add to global vars
pgv.match_intcpt = 0.0
pgv.match_factor = 0.25
pgv.max_branches = 500
pgv.max_mods = 1
pgv.modification_set = ["natural", "training_set"]
pgv.Sp_stats_plots ='n'
pgv.parallel_discovery = 'y'


start = datetime.now()
read_pytheas_file("MS_data_file")  # read in MS data to pgv.ms2_dict

pgv.n_ms2_keys = len(list(pgv.ms2_dict.keys()))
print("number of spectra = ", pgv.n_ms2_keys)

error = read_standard_files()

read_pytheas_file("rna_mod_defs_light")
build_mass_dict()
 
if __name__ == "__main__":
    print("********MAIN*********")

    with multiprocessing.Pool(initializer=init_worker, initargs=(pgv)) as pool:
    # with multiprocessing.Pool() as pool:
        ms2_file_list = pool.map(discovery_fork, [0,1,2,3,4,5])
    

# ms2_ctr = 0
# # master_match_dict = {}
# ms2_file_list = []

# if pgv.parallel_discovery == "y":
    
#     print("__name__", __name__)


# else:
#     pgv.ms2_file_list = discover_spectra()

# pgv.discovery_dict = {}
# for ms2_file in pgv.ms2_file_list:
#     ms2_key = int(ms2_file.split("/")[-1].split(".")[0].split("_")[-1])
#     pgv.discovery_dict[ms2_key] = load_pickle(ms2_file)

# pgv.unpacked_discovery_dict = unpack_master_discovery_dict(pgv.discovery_dict)

# discovery_job = pgv.discovery_job.split("_")[-1]
# discovery_file = "discovery_output" + "_" + discovery_job
# output_match_dict_file(pgv.unpacked_discovery_dict, discovery_file)

# pgv.top_discovery_dict, pgv.seq_discovery_dict, pgv.discovery_dict = consolidated_match_output(pgv.unpacked_discovery_dict, "consolidated_discovery_output")
 
# #TODO  add Sp to validation output
# # validate_discovery()

# json_dir = os.path.join(pgv.job_dir, "pytheas_json_files")
# Path(json_dir).mkdir(parents=True, exist_ok=True)
 
# save_json_files(pgc.discovery_json, json_dir)

# # TODO try to match on sequence
# #     if pgv.plot_sequence_map == 'y':
# #         max_seq_len = max([len(mdict["seq_list"]) for mdict in pgv.mol_dict.values()])
# #         if max_seq_len < 100:
# #             ma.make_sequence_plot("discovery_sequence_map")
# #         else:
# #             ma.make_long_sequence_plot("discovery_sequence_map")

# print()
# print("discovery took :", datetime.now() - start)

# Discovery()


