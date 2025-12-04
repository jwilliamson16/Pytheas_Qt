#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 14:31:53 2022

@author: jrwill
"""

# adapted from python version of isodist_3.2b   march 2017 version
import sys
# import os
# from scipy.fft import rfft, irfft,rfftfreq
# from scipy.optimize import least_squares
# import numpy as np
# from matplotlib import pyplot as plt
# import copy
# import time
from datetime import datetime

start = datetime.now()

# package_dir = "/Users/jrwill/prog/Pytheas_Folder/Pytheas_Qt_root/pytheas_dev/"
pytheas_dir =  "/Users/jrwill/prog/Pytheas_Folder/Pytheas_Qt_root/Pytheas_Qt/"
# sys.path.insert(0, package_dir)
sys.path.insert(0, pytheas_dir)


from pyteomics import mgf, mzml
import pandas as pd


from isodist_global_variables import igv
from isodist_functions import (create_atom_dict, create_model,
                               create_atom_mu, initialize_residue_mu, calc_isodist, 
                               calculate_residuals, calculate_mol_form, 
                               plot_obs_calc, find_active_residues, 
                               initialize_fit_model, fit_isotope_distribution, 
                               isodist_setup)

from minispectrum_functions import (miniSpectrum, plot_minispectrum, build_rt_indices, 
                                    build_minispec_dict, extract_minispectra, 
                                    save_rt_json, load_rt_json, extract_minispec_by_index, 
                                    fit_rt_peak, sum_spectra, quick_plot_minispectrum)

# from mod_seq_functions import parse_mod_seq, generate_mod_seq
from pytheas_global_vars import pgv,pgc
# from pytheas_objects import fragment_sequence



# STEP 1:  read MS1 datafile

step1_start = datetime.now()
MS1_data_file = "/Users/jrwill/prog/Massacre_backup/Massacre_RNA_test_set/HEK293T_28S_meth100_T1_P1-A-2_01_14293.mzML"

print("STEP 1:  reading mzML file")
igv.data = mzml.MzML(MS1_data_file)  # 30 sec

step1_end = datetime.now()
# STEP 2:  initialize minispectra from match output
print()
print("STEP 2:  initializing minispectra")
pgv.top_match_df = pd.DataFrame.from_dict(pgv.top_match_dict, orient = "index") # build dataframe from top_match_dict
igv.minispec = {idx:miniSpectrum(idx, row) for idx, row in pgv.top_match_df.iterrows()} # fast

step2_end = datetime.now()

# STEP 3:  initizalize RT indices
print()
print("STEP 3:  initializing RT indices from mzML file")
build_rt_indices(MS1_data_file) # 3 minutes -> 

#TODO   check for file, read if there, if not calculate and save
save_rt_json()
# load_rt_json()
step3_end = datetime.now()

# STEP 4:  build minispectra
print()
print("STEP 4:  extracting minispectra from mzML")
for idx, ms in igv.minispec.items():
    ms.set_mz_window()  # experimental mz array from last scan
    ms.set_rt_window()

build_minispec_dict()

#TODO  check for file, read if there, if not calculate and save
extract_minispectra()

step4_end = datetime.now()

# STEP 5:  fit RT and isotope distribution
print()
print("STEP 5:  fitting RT profiles and isotope distributions")
#TODO QC on RT peak fits   45 /1440 don't fit
for i in igv.minispec.keys():
    fit_rt_peak(i)

for i in igv.minispec.keys():
    sum_spectra(i)


# ISODIST setup

isodist_setup()

data_rows = [] # hold column data to build dataframe
ctr = 0
nspec = len(igv.minispec)
for i in igv.minispec.keys():
    print()
    print("############## fitting # ", ctr, " of ", nspec, " key = ", i)
    data_rows.append(fit_isotope_distribution(i))
    ctr += 1

isodist_df = pd.DataFrame(data_rows)
merged_df = pgv.top_match_df.join(isodist_df.set_index(pgv.top_match_df.index))

merged_df.to_excel(igv.isodist_file)

step5_end = datetime.now()

print()
print("Step 1 time = ", step1_end - start)
print("Step 2 time = ", step2_end - step1_end)
print("Step 3 time = ", step3_end - step2_end)
print("Step 4 time = ", step4_end - step3_end)
print("Step 5 time = ", step5_end - step4_end)
print("TOTAL TIME = ", datetime.now() - start)
# merged_df.to_excel(igv.isodist_file, index=False)
# igv.csvfile.close()

# total time is ~39 minutes
