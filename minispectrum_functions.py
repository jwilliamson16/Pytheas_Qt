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


class miniSpectrum:
    def __init__(self, index, match_row):
        for key, val in match_row.items():
            setattr(self, key, val)
        self.mz = match_row["m/z"]
        self.RT = match_row["RT"]  # RT in minutes
        self.rt = self.RT * 60.0   # rt in seconds
        rtpad = 60.0
        self.rtlo = self.rt - rtpad
        self.rthi = self.rt + rtpad
        self.index = index
        self.seq = match_row["Pytheas_sequence"]
        self.isotope = match_row["isotope"]
        self.z = abs(match_row["charge"])
        ctr = Counter(self.seq)
        mzpad = 2.0
        self.del_hl = 4 * (ctr["A"] + ctr["G"])/self.z
        if self.isotope == "light":
            self.mzlo = self.mz - mzpad
            self.mzhi = self.mz + self.del_hl + mzpad
        else:
            self.mzlo = self.mz - self.del_hl - mzpad
            self.mzhi = self.mz + mzpad
        self.rtlist = np.array([])
               
        
    def set_mz_window(self, mz): 
        self.mzidxlo = np.abs(mz - self.mzlo).argmin()
        self.mzidxhi = np.abs(mz - self.mzhi).argmin()
        self.mini_mz = mz[self.mzidxlo:self.mzidxhi]
        self.mz_idx = np.abs(self.mini_mz - self.mz).argmin()

    def set_rt_window(self, rt_list):
        self.rtidxlo = np.abs(rt_list - self.rtlo).argmin()
        self.rtidxhi = np.abs(rt_list - self.rthi).argmin()
        self.mini_rt = rt_list[self.rtidxlo:self.rtidxhi]
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


def read_pythid_file(file):
    df = pd.read_csv(file)
    print("length = ", len(df))
    print("columns = ", df.columns)
    return df


def gaussian(x, amplitude, mean, std_dev, offset):
        return offset + amplitude * np.exp(-((x - mean)**2) / (2 * std_dev**2))

def plot_minispectrum(ms):
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4) # 4x4 grid
    
    ax_main = fig.add_subplot(gs[1:4, 0:3]) # Takes up rows 1-3, columns 0-2
    ax_main.contour(ms.mini_mz, ms.rtlist, ms.ms, cmap = "hot_r")
    # fig.colorbar(contour, ax=ax_main, shrink=0.7) # Add colorbar
    ax_main.set_xlabel("m/z")
    ax_main.set_ylabel("RT")
    
    ax_plotx = fig.add_subplot(gs[0, 0:3], sharex=ax_main) # Top row, shares x-axis with main
    ax_plotx.plot(ms.mini_mz, ms.rt_slice)
    ax_plotx.tick_params(axis='x', labelbottom=False) # Hide x-axis labels for top hist

    ax_ploty = fig.add_subplot(gs[1:4, 3], sharey=ax_main) # Rightmost column, shares y-axis with main
    ax_ploty.plot(ms.mz_slice, ms.rtlist)
    ax_ploty.tick_params(axis='y', labelleft=False) # Hide x-axis labels for top hist

    plt.tight_layout()
    plt.show()

def initialize_minispec(file, mz, rt_list):
    match_df = read_pythid_file(file)
    minispec = {idx:miniSpectrum(idx, row) for idx, row in match_df.iterrows()} # fast
    for idx, ms in minispec.items():
        ms.set_mz_window(mz)
        ms.set_rt_window(rt_list)
    return match_df, minispec

def build_rt_indices(ms_data_file):
    start = datetime.now()
    print()
    print("Building RT indices")
    
    rt_index_dict = {}
    rt_scan_dict = {}
    rt_list = []
    
    with mzml.read(ms_data_file) as reader:
       ctr = 0
       eof = 0
       while eof == 0:
            
            try:
                scan = next(reader)
            except:
                eof = 1
                continue
            if "scanList" in scan.keys():
    
                rt = float(scan["scanList"]["scan"][0]["scan start time"])
                rt_list.append(rt)
                rt_index_dict[ctr] = rt  # index to search thru for minispectra
                rt_scan_dict[rt] = scan["id"]  # scan index to extract minispectra
            
            if ctr%100 == 0:
                print(ctr, datetime.now() - start)
            ctr += 1
    
    print("index took: ", datetime.now() - start)
    mz = scan["m/z array"]
    rt_list = np.array(rt_list)
    return rt_index_dict, rt_scan_dict, mz, rt_list


def build_mini_spec_dict(minispec, rt_index_dict):
    print()
    print("Building minispectrum dictionary")
    rt_dict = {}
    pidx_dict = {}
    
    start = datetime.now()
    for idx, rt in rt_index_dict.items():
        # RT = round(float(rt)/60.0,2)
        for pidx, ms in minispec.items():
            if ms.rtlo <= rt <= ms.rthi:
                if rt in rt_dict:
                    rt_dict[rt].append(pidx)
                else:
                    rt_dict[rt] = [pidx]
                if pidx in pidx_dict:
                    pidx_dict[pidx].append(rt)
                else:
                    pidx_dict[pidx] = [rt]
                    
                    
    print("rt_dict took: ", datetime.now()-start)  # 2.2 sec
    
    pidx_list = np.array(list(pidx_dict.keys()))

    return rt_dict, pidx_list

def extract_minispectra(data, minispec, rt_scan_dict, rt_dict, rt_list):
    print()
    print("Extracting minispectra")

    start = datetime.now()
    ctr = 0
    n = len(rt_list)
    for rt in rt_list:
        if rt not in rt_dict:
            continue
        st = datetime.now()
        scan_id = rt_scan_dict[rt]
        ctr += 1
        plist = rt_dict[rt]
        scan = data.get_by_id(scan_id)
        intensity = scan["intensity array"]
        for p in plist:
            minispec[p].add_row(intensity, rt)
        print("rt: ", rt, ctr, " of ", n, " #pidx = ", len(rt_dict[rt]), "t = ", datetime.now() - st)
    
    for idx, ms in minispec.items():
        ms.projections()

    print("extracting minispectra took: ", datetime.now() - start)


def build_minispectra(pythid_file, ms_data_file):
    start = datetime.now()  # 30 sec
    print("reading MZML file")
    data = mzml.MzML(ms_data_file)
    print(" open mzml took ", datetime.now() - start)
        
    match_df = read_pythid_file(pythid_file) # should get this from match
    minispec = {idx:miniSpectrum(idx, row) for idx, row in match_df.iterrows()} # fast
    
    rt_index_dict, rt_scan_dict, mz, rt_list = build_rt_indices(ms_data_file)
    
    for idx, ms in minispec.items():
        ms.set_mz_window(mz)
        ms.set_rt_window(rt_list)
    
    rt_dict, pidx_list = build_mini_spec_dict(minispec, rt_index_dict)
    
    extract_minispectra(data, minispec, rt_scan_dict, rt_dict, rt_list)
    
    return match_df, minispec

# pythid_file = "/Users/jrwill/prog/Pytheas_Folder/Pytheas_data/human_28S_RNA_for_massacre/Match_Job_099/HEK293T_28S_meth100_T1_P1-A-2_01_14293_massacre.csv"
# ms_data_file = "/Users/jrwill/prog/Massacre_backup/Massacre_RNA_test_set/HEK293T_28S_meth100_T1_P1-A-2_01_14293.mzML"

# match_df, minispec = build_minispectra(pythid_file, ms_data_file)





