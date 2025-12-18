#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 14:31:53 2022

@author: jrwill
"""

# adapted from python version of isodist_3.2b   march 2017 version
import sys
import os
# from scipy.fft import rfft, irfft,rfftfreq
# from scipy.optimize import least_squares
# import numpy as np
# from matplotlib import pyplot as plt
# import copy
# import time
from datetime import datetime
from pathlib import Path
from copy import deepcopy

# package_dir = "/Users/jrwill/prog/Pytheas_Folder/Pytheas_Qt_root/pytheas_dev/"
pytheas_dir =  "/Users/jrwill/prog/Pytheas_Folder/Pytheas_Qt_root/Pytheas_Qt/"
# sys.path.insert(0, package_dir)
sys.path.insert(0, pytheas_dir)


from pyteomics import mgf, mzml
import pandas as pd
import numpy as np

from pytheas_IO import read_pytheas_file, save_pickle_set, load_pickle_set

# from isodist_global_variables import igv
from isodist_functions import (create_atom_dict, create_model,
                               create_atom_mu, initialize_residue_mu, calc_isodist, 
                               calculate_residuals, calculate_mol_form, 
                               plot_obs_calc, find_active_residues, 
                               initialize_fit_model, fit_isotope_distribution, 
                               isodist_setup, Load_thru_RT, Load_thru_Minispec,
                               Load_thru_RT_fit)

from minispectrum_functions import (miniSpectrum, plot_minispectrum, build_rt_indices, 
                                    build_minispec_dict, extract_minispectra,
                                    fit_rt_peak, sum_spectra)

# from mod_seq_functions import parse_mod_seq, generate_mod_seq
from pytheas_global_vars import pgv,pgc
# from pytheas_objects import fragment_sequence


def isodist():

    # STEP 1:  read MS1 datafile
    
    start = datetime.now()

    # MS1_data_file = "/Users/jrwill/prog/Massacre_backup/Massacre_RNA_test_set/HEK293T_28S_meth100_T1_P1-A-2_01_14293.mzML"
    
    # print("STEP 1:  reading mzML file")
    # pgv.data = mzml.MzML(MS1_data_file)  # 30 sec
    # read_pytheas_file(pgv.MS1_data_file)
    
    # step1_end = datetime.now()
    # STEP 1:  initialize minispectra from match output
    print()
    print("STEP 1:  initializing minispectra")
    pgv.top_match_df = pd.DataFrame.from_dict(pgv.top_match_dict, orient = "index") # build dataframe from top_match_dict
    pgv.minispec = {idx:miniSpectrum(idx, row) for idx, row in pgv.top_match_df.iterrows()} # fast
    
    step1_end = datetime.now()
    
    # STEP 2:  initizalize RT indices
    print()
    print("STEP 2:  initializing RT indices from mzML file")
    # read file into data iterator, then iterate thru
    
    # rt_json_dir =  os.path.join(pgv.job_dir, "isodist_json_files")
    # rt_json_dir =  pgv.job_dir
    if pgv.load_rt_dict == "y":
        Load_thru_RT()
        
    else: 
        build_rt_indices() # 3 minutes -> 
        # pgv.rt_list_dict = {"rt_list": list(pgv.rt_list)}
        # pgv.mz_exp_dict = {"mz_exp": list(pgv.mz_exp)}
        save_pickle_set(pgc.isodist_rt_pickle, pgv.job_dir)

    #TODO   check for file, read if there, if not calculate and save
    
    step2_end = datetime.now()
    
    # STEP3:  build minispectra
    print()
    print("STEP 3:  extracting minispectra from mzML")
    for idx, ms in pgv.minispec.items():
        ms.set_mz_window()  # experimental mz array from last scan
        ms.set_rt_window()
    
    build_minispec_dict()
    
    #TODO  check for file, read if there, if not calculate and save
    if pgv.load_minispectra == "y":
        Load_thru_Minispec()
    else:
        extract_minispectra()
        save_pickle_set(pgc.isodist_minispec_pickle, pgv.job_dir)
      
    step3_end = datetime.now()
    
    # STEP 4:  fit RT profile
    print()
    print("STEP 4:  fitting RT profiles")
    #TODO QC on RT peak fits   45 /1440 don't fit
    
    if pgv.load_rt_fit == "y":
        Load_thru_RT_fit()
    else:
        rt_data_rows = []
        for i in pgv.minispec.keys():
            rt_data_rows.append(fit_rt_peak(i))
        
        for i in pgv.minispec.keys():
            sum_spectra(i)
        
        pgv.minispec_rt = pgv.minispec.copy()
        save_pickle_set(pgc.isodist_rt_fit_pickle, pgv.job_dir)
        step4_end = datetime.now()
    
    
    # STEP 5:  Isodist fitting
    print()
    print("STEP5: fitting of isotope distributions")
    isodist_setup()
    
    data_rows = [] # hold column data to build dataframe
    ctr = 0
    nspec = len(pgv.minispec)
    
    if "plotfit" in pgv.isodist_plot_options:
        pgv.isodist_plot_dir = os.path.join(pgv.job_dir, "isodist_plot_dir")
        Path(pgv.isodist_plot_dir).mkdir(parents=True, exist_ok=True)
    
    for i in pgv.minispec.keys():
        print()
        print("############## fitting # ", ctr, " of ", nspec, " key = ", i)
        data_rows.append(fit_isotope_distribution(i))
        ctr += 1
    
    pgv.minispec_isodist = pgv.minispec.copy()
    save_pickle_set(pgc.isodist_fit_pickle, pgv.job_dir)
    
    isodist_df = pd.DataFrame(data_rows)
    rt_df = pd.DataFrame(rt_data_rows)
    merged_df = pgv.top_match_df.join(rt_df.set_index(pgv.top_match_df.index))
    merged_df2 = merged_df.join(isodist_df.set_index(pgv.top_match_df.index))
    
    isodist_job = pgv.job_dir.split("_")[-1]
    isodist_output = os.path.join(pgv.job_dir, "isodist_output_" + isodist_job + ".xlsx")
    merged_df2.to_excel(isodist_output)
    
    step5_end = datetime.now()
    
    print()
    print("Step 1 time (initialize minispec) = ", step1_end - start)
    print("Step 2 time (initialize RT indices) = ", step2_end - step1_end)
    print("Step 3 time = (extract minispectra)", step3_end - step2_end)
    print("Step 4 time (fit RT profiles) = ", step4_end - step3_end)
    print("Step 5 time (fit isotope distributions) = ", step5_end - step4_end)
    print("TOTAL TIME = ", datetime.now() - start)
    # merged_df.to_excel(igv.isodist_file, index=False)
    # igv.csvfile.close()
    
    # total time is ~39 minutes
