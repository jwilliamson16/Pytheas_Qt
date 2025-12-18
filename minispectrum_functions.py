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
from matplotlib import rcParams

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


def gaussian(x, amplitude, mean, std_dev, offset):
        return offset + amplitude * np.exp(-((x - mean)**2) / (2 * std_dev**2))

def make_fit_table(ms):
    rt_pars = ["amp", "ctr", "wid", "b"]
    rt_par_labels = ["_".join(["rt", p, "fit"]) for p in rt_pars]
    rt_par_row_labels =  ["_".join(["rt", p]) for p in rt_pars]
    rt_par_err_labels = ["_".join(["rt", p, "err"]) for p in rt_pars]
    
    id_pars = pgv.parlabel
    id_par_labels = id_pars
    id_par_row_labels = id_pars
    id_par_err_labels = ["_".join([ p, "err"]) for p in id_pars]

    par_labels = id_par_labels + rt_par_labels
    par_err_labels = id_par_err_labels + rt_par_err_labels
    
    col_labels = ["FIT", "ERROR"]
    row_labels = id_par_row_labels + rt_par_row_labels

    data = [[getattr(ms,par), getattr(ms, perr)] for par, perr in zip(par_labels, par_err_labels)]
    data = np.round(np.asarray(data), decimals = 4)
    
    print("row_labels", row_labels)
    print("col_labels", col_labels)
    print("data", data)
    
    return row_labels, col_labels, data
    # row_labels = pars
    

def plot_minispectrum(ms):
    
    pfont = 'Arial'
    fs = 9    #fontsize for peak labels in points
    rcParams.update({'font.size': fs, 'font.family': "sans-serif", "font.sans-serif": pfont})

    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4) # 4x4 grid
    
    ax_main = fig.add_subplot(gs[1:4, 0:3]) # Takes up rows 1-3, columns 0-2   contour plot
    ax_main.contour(ms.mini_mz, ms.mini_rt, ms.ms, cmap = "hot_r")
    ax_main.set_xlabel("m/z")
    ax_main.set_ylabel("RT")
    
    ax_plotx = fig.add_subplot(gs[0, 0:3], sharex=ax_main) # Top row isotope distribution
    ax_plotx.plot(ms.mini_mz, ms.rt_slice, ".k")
    ax_plotx.plot(ms.mini_mz, ms.mz_fit, "-r")
    ax_plotx.tick_params(axis='x', labelbottom=False) # Hide x-axis labels for top hist

    ax_ploty = fig.add_subplot(gs[1:4, 3], sharey=ax_main) # Rightmost column ion chromatogram
    ax_ploty.plot(ms.mz_slice, ms.mini_rt, ".k")
    ax_ploty.plot(ms.fit, ms.mini_rt, "-r")
    ax_ploty.plot([0, max(ms.mz_slice)],[ms.rt, ms.rt],  "-b", label = "rt MS1")
    ax_ploty.tick_params(axis='y', labelleft=False) # Hide x-axis labels for top hist

    if ms.rt_ctr_fit != 0 :
        ax_ploty.plot([0,ms.rt_amp_fit],[ms.rt_ctr_fit, ms.rt_ctr_fit], "--g", label = "rt corrected")
        ax_ploty.plot([0, ms.rt_amp_fit], [ms.rt_ctr_fit + 2 * ms.rt_wid_fit, ms.rt_ctr_fit + 2 * ms.rt_wid_fit],  "--k", label = "+2 * sig")
        ax_ploty.plot([0, ms.rt_amp_fit],[ms.rt_ctr_fit - 2 * ms.rt_wid_fit, ms.rt_ctr_fit - 2 * ms.rt_wid_fit],  "--k", label = "-2 * sig")
        ax_ploty.legend()
    
    ax_text = fig.add_subplot(gs[0,3]) # Rightmost column, top row  fit parameters
    ax_text.axis('off')
    
    rows, columns, data = make_fit_table(ms)
    table = ax_text.table(cellText = data, colLabels = columns, rowLabels = rows, 
                          loc='center', cellLoc = "left")
    table.auto_set_font_size(False)
    table.set_fontsize(6)
    hgt = 0.1
    cellDict = table.get_celld()
    for i in range(-1,len(columns)): # row labels are in -1
         for j in range(0,len(data)+1):
            if j == 0 and i == -1:
                continue
            cellDict[(j,i)].set_height(hgt)
            cellDict[(j,i)].set_linewidth(0.1)
    
    #TODO fix minispec to have modseq
    plt.suptitle(" ".join([str(ms.ms2_key), generate_mod_seq(ms.frag3), str(ms.z), ms.mol ]))
    plt.tight_layout()
    
    plt.savefig("test_mini.png", dpi = 300)

# def quick_plot_minispectrum(ms):
#     fig = plt.figure(figsize=(8, 8))
#     gs = GridSpec(4, 4) # 4x4 grid
    
#     ax_main = fig.add_subplot(gs[1:4, 0:3]) # Takes up rows 1-3, columns 0-2
#     ax_main.contour(ms.mini_mz, ms.rtlist, ms.ms, cmap = "hot_r")
#     # fig.colorbar(contour, ax=ax_main, shrink=0.7) # Add colorbar
#     ax_main.set_xlabel("m/z")
#     ax_main.set_ylabel("RT")
    
#     ax_plotx = fig.add_subplot(gs[0, 0:3], sharex=ax_main) # Top row, shares x-axis with main
#     ax_plotx.plot(ms.mini_mz, ms.rt_slice, ".k")
#     # ax_plotx.plot(ms.mini_mz, ms.mz_fit, "-r")
#     ax_plotx.tick_params(axis='x', labelbottom=False) # Hide x-axis labels for top hist

#     ax_ploty = fig.add_subplot(gs[1:4, 3], sharey=ax_main) # Rightmost column, shares y-axis with main
#     ax_ploty.plot(ms.mz_slice, ms.rtlist, ".k")
#     # ax_ploty.plot(ms.fit, ms.mini_rt, "-r")
#     # ax_ploty.plot([0,ms.rt_amp_fit],[ms.rt_fit, ms.rt_fit], "--g", label = "rt corrected")
#     # ax_ploty.plot([0, ms.rt_amp_fit], [ms.rt_fit + 2 * ms.rt_width_fit, ms.rt_fit + 2 * ms.rt_width_fit],  "--k", label = "2 * sig")
#     # ax_ploty.plot([0, ms.rt_amp_fit],[ms.rt_fit - 2 * ms.rt_width_fit, ms.rt_fit - 2 * ms.rt_width_fit],  "--k", label = "2 * sig")
  
#     ax_ploty.tick_params(axis='y', labelleft=False) # Hide x-axis labels for top hist

#     plt.suptitle(" ".join([ms.ms2_key, generate_mod_seq(ms.frag3), str(ms.z), ms.mol]))
#     plt.tight_layout()
#     plt.show()

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

def extract_minispectra():
    print()
    print("Extracting minispectra")

    read_pytheas_file("MS1_data_file")

    start = datetime.now()
    ctr = 0
    n = len(pgv.rt_list)
    
    
    key_type = type(list(pgv.rt_scan_dict.keys())[0])  # determine type of rt_scan_dict_keys

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


def fit_rt_peak(idx):
    pars = ["amp", "ctr", "wid", "b"]
    par_labels = ["_".join(["rt", p, "fit"]) for p in pars]
    par_err_labels = ["_".join(["rt", p, "err"]) for p in pars]
    par_relerr_labels = ["_".join(["rt", p, "re"]) for p in pars]
    rt_x = pgv.minispec[idx].mini_rt
    rt_y = pgv.minispec[idx].mz_slice
    
    #TODOO fix this in miniSpectrum
    if len(rt_x) > len(rt_y):
        rt_x = rt_x[0:len(rt_y)]
    if len(rt_y) > len(rt_x):
        rt_y = rt_y[0:len(rt_x)]
    rt_y = rt_y[0:len(rt_x)]
    
    amp = max(rt_y)
    ctr = pgv.minispec[idx].rt
    print("initial ctr = ", ctr)
    sig = 1.0
    baseline = 0
    initial_guesses = [amp, ctr, sig, baseline]

    try:
        params, covariance = curve_fit(gaussian, rt_x, rt_y, p0=initial_guesses)
    except:
        params = np.zeros(4)
        covariance = np.zeros((4,4))
        
    if params[1] == 0: # refit 
        
        amp = max(rt_y)
        ctr_idx = np.where(pgv.minispec[idx].mz_slice == amp)[0][0]
        ctr = pgv.minispec[idx].mini_rt[ctr_idx]
        print("refitting, ctr = ", ctr)
        sig = 1.0
        baseline = 0
        initial_guesses = [amp, ctr, sig, baseline]

        try:
            params, covariance = curve_fit(gaussian, rt_x, rt_y, p0=initial_guesses)
        except:
            params = np.zeros(4)
            covariance = np.zeros((4,4))

    fit = gaussian(rt_x, *params)
    residuals = rt_y - fit
    chi_squared = np.sum((residuals)**2) #
    pgv.minispec[idx].fit = fit
    print("index", idx, "rt fit params", params)
    amp, ctr, sig, baseline = params
    
    
    perr = np.sqrt(np.diag(covariance))
    if np.inf in perr:
        perr = np.zeros(4)
    
    rel_err = np.zeros(4)
    for i in range(len(params)):
        if params[i] != 0:
            rel_err[i] = perr[i]/params[i]
            
    pars = ["amp", "ctr", "wid", "b"]
    
    for p, v in zip(par_labels, params):
        setattr(pgv.minispec[idx], p, v)
    for p, v in zip(par_err_labels, perr):
        setattr(pgv.minispec[idx], p, v)
    for p, v in zip(par_relerr_labels, rel_err):
       setattr(pgv.minispec[idx], p, v)

    column_dict = {col:dat for col, dat in zip(par_labels + ["rt_chisq"] + par_err_labels + par_relerr_labels, params.tolist() + [chi_squared] + perr.tolist() + rel_err.tolist())}
    
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
    if ms.rt_ctr_fit == 0:
        ms.rt_slice_sum = ms.rt_slice
    else:
        rt_lo = ms.rt_ctr_fit - 2 * ms.rt_wid_fit
        rt_hi = ms.rt_ctr_fit + 2 * ms.rt_wid_fit
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
    

