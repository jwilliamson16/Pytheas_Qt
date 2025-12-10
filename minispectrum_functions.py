#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 19:39:10 2025

@author: jrwill
"""

from collections import Counter
from datetime import datetime


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
from pyteomics import mgf, mzml

from pytheas_global_vars import pgv, pgc
# from isodist_global_variables import pgv
from mod_seq_functions import generate_mod_seq
from pytheas_IO import read_pytheas_file

class miniSpectrum:
    def __init__(self, index, match_row):
        for key, val in match_row.items():
            # print(key, val)
            setattr(self, str(key), val)
        self.mz1 = match_row["mz1"]
        self.rt = match_row["rt"]  # rt in seconds internal in pytheas
        self.RT = self.rt/60.0   # rt in minutes
        self.rtlo = self.rt - pgc.rtpad
        self.rthi = self.rt + pgc.rtpad
        self.index = index
        self.frag3 = match_row["frag3"]
        # self.seq = match_row["Pytheas_sequence"]
        self.label = match_row["label"]
        self.z = abs(match_row["z"])
        self.ufrag = match_row["seq_list"][0] # unique_frag_dict key
        if pgv.ion_mode == "-":
            zs = -self.z
        else:
            zs = self.z
            
        mz_light = pgv.unique_frag_dict[self.ufrag]["ion_frag_dict"]["light"]["mz1"][zs]
        mz_heavy = pgv.unique_frag_dict[self.ufrag]["ion_frag_dict"]["heavy"]["mz1"][zs]
            
        
        # self.del_hl = 4 * (ctr["A"] + ctr["G"])/self.z  # 4 units per purine 15N-shift
        self.del_hl = mz_heavy - mz_light
        if self.label == "light":
            self.mzlo = self.mz1 - pgc.mzpad
            self.mzhi = self.mz1 + self.del_hl + pgc.mzpad
        else:
            self.mzlo = self.mz1 - self.del_hl - pgc.mzpad
            self.mzhi = self.mz1 + pgc.mzpad
        self.rtlist = np.array([]) # empty
               
        
    def set_mz_window(self): 
        self.mzidxlo = np.abs(pgv.mz_exp - self.mzlo).argmin()
        self.mzidxhi = np.abs(pgv.mz_exp - self.mzhi).argmin()
        self.mini_mz = pgv.mz_exp[self.mzidxlo:self.mzidxhi]
        self.mz_idx = np.abs(self.mini_mz - self.mz1).argmin()

    def set_rt_window(self):
        self.rtidxlo = np.abs(pgv.rt_list - self.rtlo).argmin()
        self.rtidxhi = np.abs(pgv.rt_list - self.rthi).argmin()
        self.mini_rt = pgv.rt_list[self.rtidxlo:self.rtidxhi]
        self.rt_idx = np.abs(self.mini_rt - self.rt).argmin()

    def add_row(self, intensity, rt):
        try:
            self.ms = np.append(self.ms,[intensity[self.mzidxlo:self.mzidxhi]], axis = 0)
            self.rtlist = np.append(self.rtlist, rt)
        except:
            self.ms = np.array([intensity[self.mzidxlo:self.mzidxhi]])
            self.rtlist = np.array(rt)
        
    def projections(self):
        self.mz_slice = self.ms[:,self.mz_idx]
        self.rt_slice = self.ms[self.rt_idx,:]
        self.mini_rt = self.rtlist

#TODO mini_mz and rt_list are not the same???

# def read_pythid_file(file):
#     df = pd.read_csv(file)
#     print("length = ", len(df))
#     print("columns = ", df.columns)
#     return df


def gaussian(x, amplitude, mean, std_dev, offset):
        return offset + amplitude * np.exp(-((x - mean)**2) / (2 * std_dev**2))

# def save_rt_json():
#     with open("rt_index_dict.json", "w") as json_file:
#         json.dump(pgv.rt_index_dict, json_file, indent=4)
#     with open("rt_scan_dict.json", "w") as json_file:         # this may convert index to string...
#         json.dump(pgv.rt_scan_dict, json_file, indent=4)
#     rt_list_dict = {"rt_list": list(pgv.rt_list)}
#     mz_exp_dict = {"mz_exp": list(pgv.mz_exp)}
#     with open("rt_list_dict.json", "w") as json_file:
#         json.dump(rt_list_dict, json_file, indent=4)
#     with open("mz_exp_dict.json", "w") as json_file:
#         json.dump(mz_exp_dict, json_file, indent=4)

# def load_rt_json():
#     with open('rt_index_dict.json', 'r') as file:
#         pgv.rt_index_dict = json.load(file)
#     with open('rt_scan_dict.json', 'r') as file:
#         pgv.rt_scan_dict = json.load(file)
#     with open('rt_list_dict.json', 'r') as file:
#         rt_list_dict = json.load(file)
#     with open('mz_exp_dict.json', 'r') as file:
#         mz_exp_dict = json.load(file)
#         pgv.rt_list = np.array(rt_list_dict["rt_list"])
#         pgv.mz_exp = np.array(mz_exp_dict["mz_exp"])
        


def plot_minispectrum(ms):
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4) # 4x4 grid
    
    ax_main = fig.add_subplot(gs[1:4, 0:3]) # Takes up rows 1-3, columns 0-2
    ax_main.contour(ms.mini_mz, ms.rtlist, ms.ms, cmap = "hot_r")
    # fig.colorbar(contour, ax=ax_main, shrink=0.7) # Add colorbar
    ax_main.set_xlabel("m/z")
    ax_main.set_ylabel("RT")
    
    ax_plotx = fig.add_subplot(gs[0, 0:3], sharex=ax_main) # Top row, shares x-axis with main
    ax_plotx.plot(ms.mini_mz, ms.rt_slice_sum, ".k")
    ax_plotx.plot(ms.mini_mz, ms.mz_fit, "-r")
    ax_plotx.tick_params(axis='x', labelbottom=False) # Hide x-axis labels for top hist

    ax_ploty = fig.add_subplot(gs[1:4, 3], sharey=ax_main) # Rightmost column, shares y-axis with main
    ax_ploty.plot(ms.mz_slice, ms.rtlist, ".k")
    ax_ploty.plot(ms.fit, ms.mini_rt, "-r")
    ax_ploty.plot([0,ms.rt_amp_fit],[ms.rt_fit, ms.rt_fit], "--g", label = "rt corrected")
    ax_ploty.plot([0, ms.rt_amp_fit], [ms.rt_fit + 2 * ms.rt_width_fit, ms.rt_fit + 2 * ms.rt_width_fit],  "--k", label = "2 * sig")
    ax_ploty.plot([0, ms.rt_amp_fit],[ms.rt_fit - 2 * ms.rt_width_fit, ms.rt_fit - 2 * ms.rt_width_fit],  "--k", label = "2 * sig")
  
    ax_ploty.tick_params(axis='y', labelleft=False) # Hide x-axis labels for top hist

    plt.suptitle(" ".join([ms.Pytheas_sequence, str(ms.charge), ms.molecule_ID ]))
    plt.tight_layout()
    plt.show()

def quick_plot_minispectrum(ms):
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4) # 4x4 grid
    
    ax_main = fig.add_subplot(gs[1:4, 0:3]) # Takes up rows 1-3, columns 0-2
    ax_main.contour(ms.mini_mz, ms.rtlist, ms.ms, cmap = "hot_r")
    # fig.colorbar(contour, ax=ax_main, shrink=0.7) # Add colorbar
    ax_main.set_xlabel("m/z")
    ax_main.set_ylabel("RT")
    
    ax_plotx = fig.add_subplot(gs[0, 0:3], sharex=ax_main) # Top row, shares x-axis with main
    ax_plotx.plot(ms.mini_mz, ms.rt_slice, ".k")
    # ax_plotx.plot(ms.mini_mz, ms.mz_fit, "-r")
    ax_plotx.tick_params(axis='x', labelbottom=False) # Hide x-axis labels for top hist

    ax_ploty = fig.add_subplot(gs[1:4, 3], sharey=ax_main) # Rightmost column, shares y-axis with main
    ax_ploty.plot(ms.mz_slice, ms.rtlist, ".k")
    # ax_ploty.plot(ms.fit, ms.mini_rt, "-r")
    # ax_ploty.plot([0,ms.rt_amp_fit],[ms.rt_fit, ms.rt_fit], "--g", label = "rt corrected")
    # ax_ploty.plot([0, ms.rt_amp_fit], [ms.rt_fit + 2 * ms.rt_width_fit, ms.rt_fit + 2 * ms.rt_width_fit],  "--k", label = "2 * sig")
    # ax_ploty.plot([0, ms.rt_amp_fit],[ms.rt_fit - 2 * ms.rt_width_fit, ms.rt_fit - 2 * ms.rt_width_fit],  "--k", label = "2 * sig")
  
    ax_ploty.tick_params(axis='y', labelleft=False) # Hide x-axis labels for top hist

    plt.suptitle(" ".join([generate_mod_seq(ms.frag3), str(ms.z), ms.mol]))
    plt.tight_layout()
    plt.show()

# def initialize_minispec(file, mz, rt_list):
#     match_df = read_pythid_file(file)
#     minispec = {idx:miniSpectrum(idx, row) for idx, row in match_df.iterrows()} # fast
#     for idx, ms in minispec.items():
#         ms.set_mz_window(mz)
#         ms.set_rt_window(rt_list)
#     return match_df, minispec

def build_rt_indices():
    start = datetime.now()
    print()
    print("Building RT indices")
    
    read_pytheas_file("MS1_data_file") 
    pgv.rt_index_dict = {}
    pgv.rt_scan_dict = {}
    pgv.rt_list = []
    
    # with mzml.read(ms_data_file) as reader:
    ctr = 0
    for scan in pgv.MS1_data:
        
        if ctr == 0:
            pgv.mz_exp = scan["m/z array"]

        if "scanList" in scan.keys():

            rt = float(scan["scanList"]["scan"][0]["scan start time"])
            pgv.rt_list.append(rt)
            pgv.rt_index_dict[ctr] = rt  # index to search thru for minispectra
            pgv.rt_scan_dict[rt] = scan["id"]  # scan index to extract minispectra
        
        if ctr%100 == 0:
            print(ctr, datetime.now() - start)
        ctr += 1
        

    print("index took: ", datetime.now() - start)
    # pgv.mz_exp = scan["m/z array"]
    pgv.rt_list = np.array(pgv.rt_list)
    # return rt_index_dict, rt_scan_dict, mz, rt_list



       # eof = 0
       # while eof == 0:
            
       #      try:
       #          scan = next(reader)
       #      except:
       #          eof = 1
       #          continue
       #      if "scanList" in scan.keys():
    
       #          rt = float(scan["scanList"]["scan"][0]["scan start time"])
       #          pgv.rt_list.append(rt)
       #          pgv.rt_index_dict[ctr] = rt  # index to search thru for minispectra
       #          pgv.rt_scan_dict[rt] = scan["id"]  # scan index to extract minispectra
            
       #      if ctr%100 == 0:
       #          print(ctr, datetime.now() - start)
       #      ctr += 1
    
   

    


def build_minispec_dict():
    print()
    print("Building minispectrum dictionary")
    pgv.rt_dict = {}
    pidx_dict = {}
    
    start = datetime.now()
    for idx, rt in pgv.rt_index_dict.items():
        # RT = round(float(rt)/60.0,2)
        for pidx, ms in pgv.minispec.items():  # pidx is the index key, not row #
            if ms.rtlo <= rt <= ms.rthi:
                if rt in pgv.rt_dict:
                    pgv.rt_dict[rt].append(pidx)
                else:
                    pgv.rt_dict[rt] = [pidx]
                if pidx in pidx_dict:
                    pidx_dict[pidx].append(rt)
                else:
                    pidx_dict[pidx] = [rt]
                    
                    
    print("rt_dict took: ", datetime.now()-start)  # 2.2 sec
    
    pgv.pidx_list = np.array(list(pidx_dict.keys()))

    # return rt_dict, pidx_list

def extract_minispectra():
    print()
    print("Extracting minispectra")

    read_pytheas_file("MS1_data_file")
    # pgv.MS1_data = mzml.MzML(MS1_data_file)  # 30 sec
    start = datetime.now()
    ctr = 0
    n = len(pgv.rt_list)
    
    # determine type of rt_scan_dict_keys
    key_type = type(list(pgv.rt_scan_dict.keys())[0])

#TODO print out every 100 scans
    for rt in pgv.rt_list:
        if rt not in pgv.rt_dict:
            continue
        st = datetime.now()
        if key_type == float:   
            scan_id = pgv.rt_scan_dict[rt]
        else:
            scan_id = pgv.rt_scan_dict[str(rt)]
        ctr += 1
        plist = pgv.rt_dict[rt]
        scan = pgv.MS1_data.get_by_id(scan_id)
        intensity = scan["intensity array"]
        for p in plist:
            pgv.minispec[p].add_row(intensity, rt)
        print("rt: ", rt, ctr, " of ", n, " #pidx = ", len(pgv.rt_dict[rt]), "t = ", datetime.now() - st)
    
    for idx, ms in pgv.minispec.items():
        ms.projections()

    print("extracting minispectra took: ", datetime.now() - start)


# def build_minispectra(pythid_file, ms_data_file):
# def build_minispectra(ms_data_file):
#     start = datetime.now()  # 30 sec
#     print("reading MZML file")
#     pgv.data = mzml.MzML(ms_data_file)
#     print(" open mzml took ", datetime.now() - start)
        
#     # match_df = read_pythid_file(pythid_file) # should get this from match
#     pgv.top_match_df = pd.DataFrame.from_dict(pgv.top_match_dict, orient = "index")
#     pgv.minispec = {idx:miniSpectrum(idx, row) for idx, row in pgv.top_match_df.iterrows()} # fast
    
#     pgv.rt_index_dict, pgv.rt_scan_dict, pgv.mz_exp, pgv.rt_list = build_rt_indices(ms_data_file) # 3 minutes
    
#     for idx, ms in pgv.minispec.items():
#         ms.set_mz_window(pgv.mz_exp)  # experimental mz array from last scan
#         ms.set_rt_window(pgv.rt_list)
    
#     pgv.rt_dict, pgv.pidx_list = build_mini_spec_dict(pgv.minispec, pgv.rt_index_dict)
    
#     extract_minispectra(pgv.data, pgv.minispec, pgv.rt_scan_dict, pgv.rt_dict, pgv.rt_list)
    
    # return minispec

# pythid_file = "/Users/jrwill/prog/Pytheas_Folder/Pytheas_data/human_28S_RNA_for_massacre/Match_Job_099/HEK293T_28S_meth100_T1_P1-A-2_01_14293_massacre.csv"
# ms_data_file = "/Users/jrwill/prog/Massacre_backup/Massacre_RNA_test_set/HEK293T_28S_meth100_T1_P1-A-2_01_14293.mzML"

# match_df, minispec = build_minispectra(pythid_file, ms_data_file)


def fit_rt_peak(idx):
    parlabels = ["rt_amp", "rt_fit", "rt_wid", "rt_base"]
    par_err_labels = [l + "_err" for l in parlabels]
    par_relerr_labels = [l + "re" for l in parlabels]
    rt_x = pgv.minispec[idx].mini_rt
    rt_y = pgv.minispec[idx].mz_slice
    
    #TODOO fix this in miniSpectrum
    if len(rt_x) > len(rt_y):
        rt_x = rt_x[0:len(rt_y)]
    if len(rt_y) > len(rt_x):
        rt_y = rt_y[0:len(rt_x)]
    rt_y = rt_y[0:len(rt_x)]
    # print("x", len(rt_x), "y", len(rt_y))
    amp = max(rt_y)
    ctr = pgv.minispec[idx].rt
    sig = 1.0
    baseline = 0
    # gaussian(x, amplitude, mean, std_dev, offset)
    initial_guesses = [amp, ctr, sig, baseline]
    try:
        params, covariance = curve_fit(gaussian, rt_x, rt_y, p0=initial_guesses)
    except:
        params = np.zeros(4)
        covariance = np.zeros((4,4))
    fit = gaussian(rt_x, *params)
    pgv.minispec[idx].fit = fit
    print("index", idx, "rt fit params", params)
    amp, ctr, sig, baseline = params
    
    perr = np.sqrt(np.diag(covariance))
    if np.inf in perr:
        perr = np.zeros(4)
    
    rel_err = np.zeros(4)
    for i in range(len(params)):
        if perr[i] != 0:
            rel_err[i] = params[i]/perr[i]
            
    pgv.minispec[idx].rt_amp_fit, pgv.minispec[idx].rt_fit, pgv.minispec[idx].rt_width_fit, pgv.minispec[idx].baseline = params

    column_dict = {col:dat for col, dat in zip(parlabels + par_err_labels + par_relerr_labels, params.tolist() + perr.tolist() + rel_err.tolist())}
    
    return column_dict

    #TODO make plots optional

    # fig, ax = plt.subplots()
    # plt.plot(rt_x, rt_y, ".k", label = "ion chromatogram")
    # plt.plot(rt_x, fit, "-r", label = "fit")
    # plt.xlabel("Retention time (s)")
    # plt.ylabel("Intensity")
    # plt.plot([ctr, ctr], [0,amp], "--g", label = "rt corrected")
    # plt.plot([ctr - 2 * sig, ctr - 2 * sig], [0, amp], "--k", label = "2 * sig")
    # plt.plot([ctr + 2 * sig, ctr + 2 * sig], [0, amp], "--k", label = "2 * sig")
    
    # plt.title(generate_mod_seq(pgv.minispec[idx].frag3))
    # plt.legend()
    
    
def sum_spectra(idx):
    ms = pgv.minispec[idx]
    if ms.rt_fit == 0:
        ms.rt_slice_sum = ms.rt_slice
    else:
        rt_lo = ms.rt_fit - 2 * ms.rt_width_fit
        rt_hi = ms.rt_fit + 2 * ms.rt_width_fit
        rt_lo_idx = np.abs(ms.mini_rt - rt_lo).argmin()
        rt_hi_idx = np.abs(ms.mini_rt - rt_hi).argmin()
        spec_sum = np.sum(ms.ms[rt_lo_idx:rt_hi_idx + 1,:], axis = 0)
        # print(rt_lo_idx, rt_hi_idx)
        ms.rt_slice_sum = spec_sum
  
    # plt.plot(ms.rt_slice_sum)

    
# def extract_minispec_by_index(idx):
#     rt_idx_list = [rt for rt, plist in pgv.rt_dict.items() if idx in plist]
#     for rt in rt_idx_list:
#         scan_id = pgv.rt_scan_dict[str(rt)]
#         # plist = pgv.rt_dict[rt]
#         scan = pgv.data.get_by_id(scan_id)
#         intensity = scan["intensity array"]
#         pgv.minispec[idx].add_row(intensity, rt)
    
#     pgv.minispec[idx].projections()
#         # print("rt: ", rt, ctr, " of ", n, " #pidx = ", len(pgv.rt_dict[rt]), "t = ", datetime.now() - st)
    

