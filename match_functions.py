#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 08:39:30 2023

@author: jrwill
"""


from datetime import datetime
import math
import bisect
import copy
import os
import xlsxwriter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from pytheas_global_vars import pgv, pgc
from pytheas_objects import fragment_sequence, MS2_spectrum
from digest_functions import generate_molecular_graph, generate_mod_seq
from scoring_functions import (ppm_range, ppm_offset, sumI_all, sumI, n_calc, L_calc,
                                consecutive_series)
from worksheet_functions import format_worksheet_columns
from ms2_plot import ms2_plot, labeled_matrix_plot

def match_spectra():
    ms2_ctr = 0
    ms2_match_dict = {}
    for ms2_key, spec in pgv.ms2_dict.items():

        ms2_match_dict[ms2_key] = {}
        top, ms2_match_dict[ms2_key] = matching(ms2_ctr, ms2_key)  # find all precursor matches for each spectrum
        ms2_ctr += 1
    return ms2_match_dict

def matching(ms2_ctr, ms2_key):
    spec = pgv.ms2_dict[ms2_key]
    mz1 = spec.mz1
    tol = spec.mz1 * pgv.MS1_ppm/1000000.0 # this might be too large....
    z = abs(spec.z)
    if z == 0:
        print("z = 0 for ms2_key, skipping... ", ms2_key, z)
        z = 1

#TODO check for faster matching

# import bisect

#     sorted_list = [5, 8, 10, 15, 20, 25]
#     lower_bound = 7
#     upper_bound = 18

#     # Find the insertion point for the lower bound
#     left_index = bisect.bisect_left(sorted_list, lower_bound)

#     # Find the insertion point for the upper bound
#     right_index = bisect.bisect_right(sorted_list, upper_bound)

#     elements_in_range = sorted_list[left_index:right_index]
#     print(elements_in_range)

    precursor_list = [pkey for pkey in pgv.unique_precursor_dict.keys() 
                 for t in pgv.precursor_isotopologue_list 
                 if abs(pgv.unique_precursor_dict[pkey]["mz1"] + float(t) * pgc.neutron_mass/z - mz1) < tol 
                 and abs(pgv.unique_precursor_dict[pkey]["z"]) == z]      
    top, fdict = match_score_rank(ms2_key, precursor_list)  # should be fragment_list
    if len(top) > 0:
        top_match = top[0]
        top_score = round(fdict[0]["Sp"],3)
    else:
        top_match = "none"
        top_score = 0.0
    print("matching:", ms2_ctr, "of", pgv.n_ms2_keys, "for mz1 =", round(spec.mz1,3), "best match of",len(top), "precursors is Sp = ", top_score, top_match)
    
    return top, fdict

def match_score_rank(ms2_key, updkey_list): # unique_precursor_dict...list of keys that match precursor
    spec = pgv.ms2_dict[ms2_key]
    prec_dict = {}
    pidx = 0
    for key in updkey_list:
        prec_dict[pidx] = copy.deepcopy(pgv.unique_precursor_dict[key])
        if "CID_ions" not in pgv.unique_precursor_dict[key]: # don't re-calculate ms2 ions  ( ~ 50% extra time if you do)
            m0 = pgv.unique_precursor_dict[key]["m0"]
            frag3 = pgv.unique_precursor_dict[key]["frag3"]
            label = pgv.unique_precursor_dict[key]["label"]
            ms2_ions = build_ion_series(m0, frag3, label)
            pgv.unique_precursor_dict[key]["CID_ions"] = ms2_ions  # need to sort
  
        prec_dict[pidx]["CID_ions"] =  pgv.unique_precursor_dict[key]["CID_ions"]
        pidx += 1
    
    matching_ms2(ms2_key, spec, prec_dict) # slowest step 4.8 sec...updates frag_dict
    scoring(ms2_key, prec_dict, pgv.CID_series)  # add Sp and  scores to ms2_match_dict
    
    prec_dict = rank_matches(prec_dict, pgv.ntop) # all matches
    top_keys = list(prec_dict.keys())
    top_sequences = [fragment_sequence(prec_dict[key]["frag3"]).frag for key in top_keys]

    return top_sequences, prec_dict

def build_ion_series(m0, frag3, label):
     G = generate_molecular_graph(frag3, label) # did this during digest... store in dict??
     ms2_ions = generate_ms2_ions(G, m0, label, frag3)     
     ms2_ions_sorted = dict(sorted(ms2_ions.items(), key=lambda item: item[1]["mz2"]))
     return ms2_ions_sorted


def next_node(G, mtot, node, z2_list, m0, length): # not needed, but keep to incorporate into digest/match
    nd = G.nodes[node] # node dict 
    base = nd["base"]
    grp = nd["group"]
    mgrp = pgv.nt_fragment_dict["light"][base].mass_dict[grp]
    mtot += mgrp
    cid_ions = []
    fl = nd["fragment_left"]
    fr = nd["fragment_right"]
    ml = mtot + nd["H_corr_left"] * pgv.hmass
    mr = m0 - ml + nd["H_corr_right"] * pgv.hmass

#   left fragment
    if fl != "none":
        if nd["main_side"] == "S": # side-chain does not split backbone
            left = nd["resno"]
        else:
            left =  nd["resno"]+ nd["fragment_left_offset"]
        for zs in z2_list:
            mz2 = (ml + zs * pgv.hmass)/abs(zs)
            if pgv.MS2_mzlow < mz2 < pgv.MS2_mzhigh: 
                cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": fl, "index": left, "z": zs})
        
#   right_fragment
    if fr != "none":
        if nd["main_side"] == "S": # side-chain does not split backbone
            right = nd["resno"]
        else:
            right = length - nd["resno"] + nd["fragment_right_offset"]  # right ions are n-resno
        for zs in z2_list:
            mz2 = (mr + zs * pgv.hmass)/abs(zs)
            if pgv.MS2_mzlow < mz2 < pgv.MS2_mzhigh: 
                cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": fr, "index": right, "z": zs})

#TODO make sure y-P and z-P are not redundant with losses
                if fr == 'y':  # add y-P loss
                    mz2 = (mr - pgv.losses_dict['y-P']["loss_mass"] + zs * pgv.hmass)/abs(zs)
                    cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": 'y-P', "index": right, "z": zs})

                if fr == 'z':  # add z-P loss
                    mz2 = (mr - pgv.losses_dict['z-P']["loss_mass"] + zs * pgv.hmass)/abs(zs)
                    cid_ions.append({"mz2": mz2, "intensity": 1.0, "series": 'z-P', "index": right, "z": zs})

    return mtot, cid_ions


def generate_ms2_ions(G, m0, label, seq3n):
    # ion_ctr = 0
    seq_len = len(seq3n) - 2  # of bases in sequence (omit end5 and end3) 
    seq3 = seq3n[1:-1]
    ms2_ions = {} # charged fragments
    # ion_key = 0
    if pgv.ion_mode == "-":  
        zsign = -1
    else:
        zsign = 1

    z2_list = [z * zsign for z in pgv.MS_charge_dict["MS2_charges"][seq_len]]

    mtot = 0.0
    ion_idx = 0
    for node in G.nodes:
        mtot, cid_ions = next_node(G, mtot, node, z2_list, m0, seq_len)
        for ion in cid_ions:
            ms2_ions[ion_idx] = ion
            ion_idx +=1
    
    # add losses
    losses = []
    for key, ion in ms2_ions.items():
         for loss, ldict in pgv.losses_dict.items():
             if ldict["parent"] == ion["series"]:
                 if ion["index"] > len(seq3):
                     continue
                 if ldict["loss_mass"] == 0:
                     base = seq3[int(ion["index"]) - 1]
                     lmass = ion["mz2"] - pgv.nt_fragment_dict[label][base].mass_dict["B"]/abs(int(ion["z"]))
                 else:
                     lmass = ion["mz2"] - ldict["loss_mass"]/abs(int(ion["z"]))
                
                 if pgv.MS2_mzlow < lmass < pgv.MS2_mzhigh: 
                     losses.append( {"mz2": lmass, "intensity": 1.0, "series": loss, "index": ion["index"], "z": ion["z"]})
                   
    for loss in losses:
        ms2_ions[ion_idx] = loss
        ion_idx += 1
    return ms2_ions

def matching_ms2(ms2_key, spec, prec_dict):

    ms2ppmo = float(pgv.MS2_ppm_offset)  
    ms2ppm = float(pgv.MS2_ppm)
 
    prec_charge = spec.z
    if prec_charge == 0:
        prec_charge = 1

    for key,precursor in prec_dict.items():
        ms2_vals_sorted = [v["mz2"] for v in prec_dict[key]["CID_ions"].values()]  
        ms2_keys_sorted = list(prec_dict[key]["CID_ions"].keys())
        n_ms2_keys = len(ms2_vals_sorted)

        match_dict = {}
        midx = 0
        for mz2, intensity in spec.ms2.items():
            ppmo = ppm_range(mz2, ms2ppmo)
            ppm = ppm_range(mz2 + ppmo , ms2ppm) 
            lower = mz2 + ppmo - ppm
            upper = mz2 + ppmo + ppm
            idxlo = max(bisect.bisect_left(ms2_vals_sorted,lower) - 2,0)
            idxhi = min( bisect.bisect_left(ms2_vals_sorted,upper) + 2,n_ms2_keys)
            
            for idx in range(idxlo,idxhi):  # check for ppm tolerance
                match_mass = ms2_vals_sorted[idx]
                calc_offset = ppm_offset(mz2, match_mass)
                if abs(calc_offset) < ms2ppm:
                    cid_idx = ms2_keys_sorted[idx]
                    cid_ion = precursor["CID_ions"][cid_idx]
                    match_dict[midx] = {"theo_mz2": cid_ion["mz2"],"theo_int": cid_ion["intensity"], "series": cid_ion["series"], 
                                        "index": cid_ion["index"], "z": cid_ion["z"], "obs_mz2": mz2, "obs_int": intensity, 
                                        "ppmo": calc_offset}
                    midx += 1
        rename_matched_ions(precursor, match_dict)

        precursor["matched_ions"] = match_dict
        precursor["ms2_key"] = ms2_key
        precursor["rt"] = pgv.ms2_dict[ms2_key].rt
        precursor["mz_exp"] = pgv.ms2_dict[ms2_key].mz1
        
        offsets = [(precursor["mz_exp"] + float(t) * pgc.neutron_mass/precursor["z"]  - precursor["mz1"])/precursor["mz_exp"] for t in pgv.precursor_isotopologue_list]
        abs_offset = [abs(o) for o in offsets]
        min_idx =abs_offset.index(min(abs_offset))
        
        precursor["offset"] = 1000000.0* offsets[min_idx]
        precursor["topolog"] = min_idx
        
def rename_matched_ions(precursor, match_dict):
    f3 = fragment_sequence(precursor["frag3"])
    seq = f3.modseq
    series_key_list = []
    key_list = []
    
    for mkey, m in match_dict.items():
        if m["series"] == "B" or m["series"] == "M-P-B" or m["series"] == "M-B":
            base = seq[m["index"] - 1].replace("[","").replace("]","")
            if m["index"] <= len(seq):
                m["series"] = m["series"].replace("B", base)
                m["index"] = ""
        series_key = m["series"] + str(m["index"]) + str(m["z"])
        if series_key not in series_key_list:
            series_key_list.append(series_key)
        else:
            key_list.append(mkey)

    for key in key_list:
        match_dict.pop(key)
        
def scoring(ms2_key, prec_dict, ion_series): # Add the Sp score to the precursor ion matches
#TODO  got rid of "all" option  fix!
    if ion_series == ['all']:
         ion_series_list = pgv.all_CID_series
    else:
        ion_series_list = ion_series
    thresh = float(pgv.MS2_peak_threshold)
    sumi_all, max_int = sumI_all(pgv.ms2_dict[ms2_key], thresh) # should this be the theo digest???? NO!

    for key, precursor in prec_dict.items():
        precursor["sumi"] = sumI(precursor["matched_ions"])
        precursor["sumi_all"] = sumi_all
        precursor["n"] = n_calc(precursor["matched_ions"], [], thresh)
        precursor["L"] = L_calc(precursor["CID_ions"], [], thresh)
        precursor["beta_list"] = [round(consecutive_series(precursor["matched_ions"],s),3) for s in ion_series_list]
        precursor["beta"] = 1.0 + sum(precursor["beta_list"])
#TODO put alpha in this eqn
        precursor["alpha"] = pgv.alpha
        if precursor["sumi_all"] == 0 or  precursor["L"] == 0:
            calc_score = 0
        else:
            calc_score = precursor["sumi"] * precursor["n"] * precursor["beta"]/(precursor["sumi_all"] * precursor["L"])
        precursor["Sp"] = calc_score
        precursor["score_details"] = [":".join([sk, str(precursor[sk])]) for sk in pgv.score_keys]

        precursor["max_int"] = max_int

        # calculate Xcorr score
        tspec = calc_theo_spectrum(precursor)
        corr_fft = np.multiply(tspec, pgv.ms2_dict[ms2_key].ft)
        corr_ift = np.fft.ifft(corr_fft)
        corr = np.real(corr_ift)
        corr = np.fft.fftshift(corr)
        xcorr = max(corr)/np.mean(corr[pgc.np2 - pgv.xcorr_avg_width: pgc.np2 - pgv.xcorr_excl_width] + corr[pgc.np2 + pgv.xcorr_excl_width: pgc.np2 + pgv.xcorr_avg_width])
        precursor["Xc"] = xcorr
        precursor["Xc_sc"] = xcorr * pgv.xcorr_length_scale * precursor["length"]
        
def calc_theo_spectrum(precursor):
        mz1 = 0.0 # is mz1 needed?
        z = precursor["z"]
        rt = 0.0
        mzarray = [mz["mz2"] for ion,mz in precursor["CID_ions"].items()]
        iarray = [ion_weight(mz["series"]) for ion,mz in precursor["CID_ions"].items()]

        tspec = MS2_spectrum( mz1, rt, z, mzarray, iarray)
        tspec.normalize_intensities(pgv.normalize_MS2)
        tspec.bin_spectrum(4096)
        tspec.fft_spectrum()
        
        return tspec.ft

def ion_weight(series):
    if "M" in series:
        series = "M"
    if series in pgc.iw_dict:
        iw = pgc.iw_dict[series]
    else:
        iw = 0.0
    return iw

def rank_matches(prec_dict, top_n):
    
    # Sp ranking
    
    ms2_matches = [s for s in list(prec_dict.keys())] # list of digest_keys
    nmatches = min(len(ms2_matches), top_n)
    nmatches = len(ms2_matches)
    if len(ms2_matches) == 0:
        return {}

    Sp_matches_sort = sorted(ms2_matches, key = lambda x: prec_dict[x]["Sp"], reverse = True)
    Sp = prec_dict[Sp_matches_sort[0]]["Sp"]

    rank = 1
    for m in Sp_matches_sort[0:nmatches]:  # only look for top hits in library
        prec_dict[m]["Sp_rank"] = rank
        if pgv.normalize_scores == 'y':
            prec_dict[m]["Sp"] = prec_dict[m]["Sp"]/Sp
            prec_dict[m]["dSp"] = 1 - prec_dict[m]["Sp"]
        else:
            prec_dict[m]["dSp"] = Sp - prec_dict[m]["Sp"]
        prec_dict[m]["dSp2"] = 0.0
        rank +=1

    if len(Sp_matches_sort) > 1:
       prec_dict[Sp_matches_sort[0]]["dSp2"] = Sp - prec_dict[Sp_matches_sort[1]]["Sp"] # difference to second best match
  
    # Xc ranking

    Xc_matches_sort = sorted(ms2_matches, key = lambda x: prec_dict[x]["Xc"], reverse = True)
    Xc = prec_dict[Xc_matches_sort[0]]["Xc"]
    
    rank = 1
    for m in Xc_matches_sort[0:nmatches]:  # only look for top hits in library
        prec_dict[m]["Xc_rank"] = rank
        if pgv.normalize_scores == 'y':
            prec_dict[m]["Xc"] =  prec_dict[m]["Xc"]/Xc
            prec_dict[m]["dXc"] = 1 - prec_dict[m]["Xc"]
        else:
            prec_dict[m]["dXc"] = Xc - prec_dict[m]["Xc"]
 
        prec_dict[m]["dXc2"] = 0.0
        rank +=1
        
    if len(Xc_matches_sort) > 1:
        prec_dict[Xc_matches_sort[0]]["dXc2"] = Xc - prec_dict[Xc_matches_sort[1]]["Xc"]
        
    if pgv.ranking == "Sp":
        sort_order =  Sp_matches_sort
    elif pgv.ranking == "Xc":
        sort_order = Xc_matches_sort
    
    ordered_dict = {idx:prec_dict[sort_order[idx]] for idx in range(len(sort_order))}
    
    return ordered_dict


def output_match_dict_file(output_file): 
              
    formatted_dict = format_match_dict()

    flist = []
    for ukey, udict in formatted_dict.items():
        if udict['Sp'] < pgv.Sp_cutoff:
            flist.append(ukey)
    
    for fkey in flist:
        formatted_dict.pop(fkey)

    udf = pd.DataFrame.from_dict(formatted_dict, orient = 'index')  # convert flat dict to dataframe for convenient excel output
    udf.fillna(99999) 
    
    fd = pgv.output_format_dict["match_output"]  # fd = "format_dict" is formatting information
    # reorder columns according to fd
    column_order = [col for col in sorted(list(fd.keys()), key = lambda x: fd[x]["order"]) if col in udf.columns] # sort order
    column_order = [col for col in column_order if fd[col]["order"] < 100] # eliminate those with 999 as "not shown"
    udf = udf[column_order]
        
    # Split the matched ions so they appear in individual columns for readability
    if "matched_ions" in udf.columns:
        udf = udf.join(pd.DataFrame(udf.pop('matched_ions').str.split(" ").tolist(), index=udf.index)
            .fillna('').rename(columns=lambda c:  "matched_ions_" + str(c + 1)))
     
    match_output_file = os.path.join(pgv.job_dir, output_file + ".xlsx")
    workbook = xlsxwriter.Workbook(match_output_file,{"nan_inf_to_errors": True})
    worksheet = workbook.add_worksheet(output_file)    
    
    for okey,odict in fd.items(): # set up  column formats for xlsxwriter, list of two for green/white alternation
        if "0" in odict["format"]:
            odict["xformat"] =[ workbook.add_format({'num_format': odict["format"], "bg_color": pgc.white}),
                                workbook.add_format({'num_format': odict["format"], "bg_color": pgc.green })]
        else:
            odict["xformat"] = [workbook.add_format({'num_format': "0", "bg_color": pgc.white}),
                                workbook.add_format({'num_format': "0", "bg_color": pgc.green})]
    bold = workbook.add_format()
    bold.set_bold()
    
    format_worksheet_columns(worksheet, udf, fd)   # sets the column widths

    for col in udf.columns:  # output header
        if "matched_" in col:
            worksheet.write(0, udf.columns.get_loc(col), col, bold)
        else:
            worksheet.write(0, udf.columns.get_loc(col), fd[col]["label"], bold)  # fd has user-defined label column if desired
    
    # write out each cell
    rowidx = 1
    keyidx = 1
    lastrow = 999
    for row in udf.iterrows():
        thisrow = row[1]["ms2_key"]
        if lastrow != thisrow:
            keyidx += 1  # keeps track of matches to the same ms2_key
        lastrow = thisrow
        colidx = 0
        for col in row[1].keys():
            if "matched_" not in col:
                worksheet.write(rowidx, colidx, row[1][col], fd[col]["xformat"][keyidx%2]) # modulo 2 for green/white
            else:
                worksheet.write(rowidx, colidx, row[1][col], fd["matched_ions"]["xformat"][keyidx%2])
            colidx +=1
        rowidx +=1
    
    workbook.close()
    
def unpack_match_dict(ms2_match_dict):
    unpacked_dict = {}
    ukey = 0
    for ms2_key, mdict in ms2_match_dict.items(): 
        for midx, match_dict in mdict.items():
            unpacked_dict[ukey] = {"ms2_key": ms2_key, "match_index": midx}
            unpacked_dict[ukey].update(match_dict)
            ukey += 1
            
    top_match_dict = {}
    if pgv.ranking == "Sp":
        rank = "Sp_rank"
    else:
        rank = "Xc_rank"

    for ukey, udict in unpacked_dict.items():
        if pgv.ranking in udict:
            if udict[rank] == 1:
                top_match_dict[ukey] = udict
 
    return unpacked_dict, top_match_dict

def format_match_dict():   # entry from the unpacked dict...reformat for table output

    formatted_dict = {}
    for midx, mdict in pgv.unpacked_match_dict.items(): # reformat lists and dicts in

        fdict = copy.deepcopy(mdict)
        f3 = fragment_sequence(fdict["frag3"]) # add various sequence representations
        fdict.update(f3.__dict__)

        for key, value in fdict.items():  # convert lists to string
            if type(value) == list:
                fdict[key] = " ".join([str(v) for v in value])
                                 
        if "matched_ions" in fdict: # convert dict to string
            fdict["matched_ions"] = reformat_matched_ions(fdict["matched_ions"])
            
        if "CID_ions" in mdict:
            fdict["CID_ions"] = reformat_theo_ions(fdict["CID_ions"])
        
        formatted_dict[midx] = fdict
        
    return formatted_dict
           
def reformat_matched_ions(ion_dict):
    ion_list = []
    for ms2_ion, idict in ion_dict.items():
        ser, idx, z = idict["series"], idict["index"], idict["z"]
        iname = str(ser) + str(idx) + "(" + str(z) + ")"
        mz_t = round(idict["theo_mz2"],3)
        mz_o = round(idict["obs_mz2"],3)
        ms2_off = round(1000000.0*(mz_t - mz_o)/mz_t, 3) # ppm
        n_int =  str(round(idict["obs_int"],2))
        ion_str = "[" + iname + "](o:" + str(mz_o) + ")(" + str(ms2_off) + "ppm)(i:" + n_int + ")(t:" + str(mz_t) + ")"
        ion_list.append(ion_str)
        
    ion_string = " ".join(ion_list)
    return ion_string

def reformat_theo_ions(ion_dict):
    ion_list = []
    for ms2_ion, idict in ion_dict.items():
        ser, idx, z = idict["series"], idict["index"], idict["z"]
        iname = str(ser) + str(idx) + "(" + str(z) + ")"
        mz_t = round(idict["mz2"],3)
        # mz_o = round(idict["obs_mz2"],3)
        # ms2_off = round(1000000.0*(mz_t - mz_o)/mz_t, 3) # ppm
        n_int =  str(round(idict["intensity"],2))
        ion_str = "[" + iname + "](mz2:" + str(mz_t) + ")(i:" + n_int + ")"
        ion_list.append(ion_str)
        
    ion_string = " ".join(ion_list)
    return ion_string

def match_output_keys(udict):
    for ukey in udict.keys():
        if len(list(udict[ukey].keys())) != 0:
            pgv.match_dict_keys = list(udict[ukey].keys())
            break

def consolidated_match_output(output_file):

    rank = pgv.ranking + "_rank"
    
    top_matches = []
    for ukey, udict in pgv.unpacked_match_dict.items():
        if rank in udict:
            if udict[rank] == 1:
                top_matches.append(ukey)
    
    print("total number of matches ", len(pgv.unpacked_match_dict.keys()))
    print("number of top matches ", len(top_matches))        
    
    dSp2_threshold = 0.1
    
    top_match_dict, seq_match_dict, seq_count_dict = {}, {}, {}
 
    # flag cases of poor discrimination

    for ukey in top_matches:
        udict = pgv.unpacked_match_dict[ukey]
        top_match_dict[ukey] = udict
        ms2_key = udict["ms2_key"]
        udict["flags"] = []
        nmatches = 0
        for key, val in pgv.unpacked_match_dict.items():
            if val["ms2_key"] == ms2_key:
                nmatches += 1
 
        if nmatches == 1:
            udict["flags"].append("unique")
        else:
            if udict["dSp2"] < dSp2_threshold:
                udict["flags"].append("dSp")
        
        frag3 = udict["frag3"]
        frag = fragment_sequence(frag3).frag
        if frag not in seq_match_dict.keys():
            ckey = 0
            seq_match_dict[frag] = {}
            seq_match_dict[frag][ckey] = udict
            seq_count_dict[frag] = 1
        else:
            ckey = seq_count_dict[frag]
            seq_match_dict[frag][ckey] = udict
            seq_count_dict[frag] += 1
        
        seq_match_dict[frag][ckey]["ukey"] = ukey  # link back to master flat match dict
 
    match_dict = {}
    
    for skey, sdict in seq_match_dict.items():
        match_dict[skey] = {"frag": skey}
        mkeys_sorted = sorted(list(sdict.keys()), key = lambda x: sdict[x]["Sp"], reverse = True)
        maxidx = mkeys_sorted[0]
     
        for key, val in sdict[maxidx].items():
            match_dict[skey][key] = val
         
        mol_list = sdict[maxidx]["seq_list"]
        mol_fr_list = [":".join(s.split(":")[0:2]) for s in mol_list]
        mol_fr_list = " ".join(mol_fr_list)
        sidx_list, rt_list, z_list, Sp_list, dSp_list = [], [], [], [], []
        for mkey in mkeys_sorted:
            mdict = sdict[mkey]
            sidx_list.append(mdict["ms2_key"])
            rt_list.append(mdict["rt"])
            z_list.append(mdict["z"])
            Sp_list.append(mdict["Sp"])
            dSp_list.append(mdict["dSp2"])
            
        match_dict[skey]["flags"] = ",".join(match_dict[skey]["flags"])
        match_dict[skey]["ms2_key_list"] = ", ".join([str(s) for s in sidx_list])
        match_dict[skey]["rt_list"] = ", ".join([str(round(r,1)) for r in rt_list])
        match_dict[skey]["z_list"] = ", ".join([str(z) for z in z_list])
        match_dict[skey]["Sp_list"] = ", ".join([str(round(s,3)) for s in Sp_list])
        match_dict[skey]["dSp_list"] = ", ".join([str(round(s,3)) for s in dSp_list])
        match_dict[skey]["n_matches"] = str(len(sidx_list))
        match_dict[skey]["mol_list"] = mol_fr_list

    fd = pgv.output_format_dict["consolidated_match"]  # formatting information
        
    cm_df = pd.DataFrame.from_dict(match_dict, orient = "index")
    column_order = [col for col in sorted(list(fd.keys()), key = lambda x: fd[x]["order"]) if col in cm_df.columns] # sort order
    column_order = [col for col in column_order if fd[col]["order"] < 100] # eliminate those with 999 as "not shown"
    cm_df = cm_df[column_order]
 
    workbook = xlsxwriter.Workbook(os.path.join(pgv.job_dir, output_file + ".xlsx"),{"nan_inf_to_errors": True})
    worksheet = workbook.add_worksheet(output_file)    
    
    for okey,odict in fd.items(): # set up  column formats for xlsxwriter
        # print("cm odict", odict)
        if "0" in odict["format"]:
            odict["xformat"] = workbook.add_format({'num_format': odict["format"]})
        else:
            odict["xformat"] = workbook.add_format({'num_format': "0"})
    bold = workbook.add_format()
    bold.set_bold()
  
    format_worksheet_columns(worksheet, cm_df, fd)   

    for col in cm_df.columns:  # output header
        if col in fd:
            if "matched_" in col:
                worksheet.write(0, cm_df.columns.get_loc(col), col, bold)
            else:
                worksheet.write(0, cm_df.columns.get_loc(col), fd[col]["label"], bold)  # fd has user-defined label column if desired
    # write out each cell
    rowidx = 1
    for row in cm_df.iterrows():
        colidx = 0
        for col in row[1].keys():
            if col in fd:
                if type(row[1][col]) == float:
                    if math.isnan(row[1][col]):
                        print("NAN", rowidx, col)
                    
                if "matched_" not in col:
                    worksheet.write(rowidx, colidx, row[1][col], fd[col]["xformat"]) 
                else:
                    worksheet.write(rowidx, colidx, row[1][col], fd["matched_ions"]["xformat"])
                colidx +=1
        rowidx +=1
        
    workbook.close()
 
    return top_match_dict, seq_match_dict, match_dict

def matrix_plot(color_matrix, row_labels, col_labels, text_dict, output_file):
    
    fs = 36  # fontsize for labels
    pfont = 'Arial'
    plt.rcParams.update({'font.size': fs, 'font.family': "sans-serif", "font.sans-serif": pfont})
    plt.rcParams.update({'font.size': fs, 'font.family': "monospace", "font.sans-serif": pfont})
    
    fig, ax = plt.subplots(layout="constrained")  # plot internal coords are 1 unit/row-col

    nr, nc, _ = color_matrix.shape
    xsize = nc + 2*max([len(s) for s in row_labels]) *fs/72 + 1
    ysize = nr + 2*max([len(s) for s in col_labels]) *fs/72 + 1
    xscale = xsize
    yscale = ysize
    
    if pgv.scale_matrix_plots == 'y':
        xtarget = 8
        ytarget = 11
        xscale = xsize/xtarget
        yscale = ysize/ytarget
        xyscale = max(xscale, yscale, 1.0)
    else:
        xyscale = 1.0
    
    print("plot size: ", xsize/xyscale, ysize/xyscale, "scale", xyscale)
    fig.set_size_inches(xsize/xyscale, ysize/xyscale) # absolute plot size for display in inches
    plt.rcParams["figure.autolayout"] = True
    plt.axis('off')
    ax.set_aspect('equal')

    labeled_matrix_plot(color_matrix, row_labels, col_labels, text_dict, fs, pgv.match_sequence_labels, xyscale, ax)
   
    fig_width, fig_height = plt.gcf().get_size_inches()
    pdffile = os.path.join(pgv.job_dir, output_file + ".pdf")
    fig.savefig(pdffile, format='pdf', dpi=300)

    plt.close()

def make_mod_text_dict(row_labels, col_labels, find_row, molecule, nc):
    mod_text_dict = {} # text for modified positions
    idx = 0 
    for mod, mdict in pgv.mod_dict.items():
        mpos = str(mdict["Position"])
        mol = mdict["Molecule"]
        if find_row:
            col_idx = col_labels.index(mpos)
            row_idx = row_labels.index(mol)
        else:
            if mol == molecule:
                midx = mdict["Position"]
                row_idx = math.floor((midx-1)/nc)
                col_idx = midx - 100 * row_idx - 1
            else:
                continue

        mbase = pgv.nt_def_dict[mdict["ID"]]["Pytheas_ID"]
        mod_text_dict[idx] = {"row": row_idx, "col": col_idx, "text": mbase}
        idx += 1
    return mod_text_dict
    
def sequence_color(n_matches, seq):
    col = pgc.dark_gray
    if n_matches == 1 and seq in pgc.natural:
        col = pgc.dark_blue
    if n_matches > 1 and seq in pgc.natural:
        col = pgc.light_blue
    if n_matches == 1 and seq not in pgc.natural:
        col = pgc.dark_green
    if n_matches > 1 and seq not in pgc.natural:
        col = pgc.light_green
    return col

def make_sequence_plot(output_file):
    
    # determine matrix dimensions nr = n_seq, nc = max_seq_len
    max_seq_len = 0
    for mol, mdict in pgv.mol_dict.items():
        slen = len(mdict["seq3"])
        if slen > max_seq_len:
            max_seq_len = slen
    n_seq = len(list(pgv.mol_dict.keys()))

    row_labels = [mol for mol in pgv.mol_dict.keys()]
    col_labels = [str(i+1) for i in range(max_seq_len)]
   
    color_matrix = np.asarray([[pgc.white_hsv for i in range(max_seq_len)] for j in range(n_seq)])
    nr, nc, _ = color_matrix.shape

    for mol, mdict in pgv.mol_dict.items(): # gray out entries with no sequence
        slen = len(mdict["seq3"])
        row_idx = row_labels.index(mol)
        for i in range(slen, max_seq_len): # gray
            color_matrix[row_idx,i] = pgc.dark_gray
    
    # text labels for mods
    mod_text_dict = make_mod_text_dict(row_labels, col_labels, True, "", 100)

    # add modifications to sequence matrix as red...to be overwritten if matched
    for key, mdict in pgv.mod_dict.items():
        molecule = mdict["Molecule"]
        idx = mdict["Position"]
        if molecule in row_labels:
            row_idx = row_labels.index(molecule)
            col_idx = idx - 1
            color_matrix[row_idx, col_idx] = pgc.red

#TODO add light color for non-top match
    for m, mdict in pgv.match_dict.items():
        mol_list = mdict["mol_list"].split(" ")
        n_matches = len(mol_list)
        seq3 = mdict["frag3"][1:-1]
        for mol_str in mol_list:
            mol, r = mol_str.split(":")
            if mol not in pgv.mol_dict:
                continue
            row_idx = row_labels.index(mol)
            
            col_fr, col_to = map(int,r.split("_"))
            seq3_idx = 0
            for col_idx in range(col_fr - 1, col_to):
                # print(seq3, len(seq3) + 1, col_idx, seq3_idx)
                if seq3_idx >= len(seq3):
                    # print("fragment too long", seq3, mol_str)
                    break
                color_matrix[row_idx, col_idx] = sequence_color(n_matches, seq3[seq3_idx])
                seq3_idx += 1

    matrix_plot(color_matrix, row_labels, col_labels, mod_text_dict, output_file)
    
def make_long_sequence_plot(output_file):
    
    # plot long sequences as blocks of 100, one plot for each sequence
    nc = 100
    for molecule, mdict in pgv.mol_dict.items():
        slen = len(mdict["seq3"])
        nr = math.floor(slen/nc)
        if slen%100 != 0:
            nr += 1
            
        row_labels = [molecule + "(" + str(100*row + 1) + ":" + str(100*(row + 1)) + ")"  for row in range(nr)]
        col_labels = [str(i+1) for i in range(nc)]
       
        color_matrix = np.asarray([[pgc.white_hsv for i in range(nc)] for j in range(nr)])
         
        for row in range(nr): # gray out cells > length of seq
            for col in range(nc):
                idx = 100*row + col + 1
                if idx > slen:
                    color_matrix[row, col] = pgc.dark_gray # dark gray
         
        # add modifications to sequence matrix...to be overwritten if matched
        for key, mdict in pgv.mod_dict.items():
            if mdict["Molecule"] != molecule:
                continue
            idx = mdict["Position"]
            row_idx = math.floor((idx-1)/nc)
            col_idx = idx - 100 * row_idx - 1
            color_matrix[row_idx,col_idx] = pgc.red

        # text labels for mods
        mod_text_dict = make_mod_text_dict(row_labels, col_labels, False, molecule, nc)
        
        # fill out sequence matrix based on unique/multiple IDs
        for m, mdict in pgv.match_dict.items():
            mol_list = mdict["mol_list"].split(" ")
            n_matches = len(mol_list)
            seq3 = mdict["frag3"][1:-1]
            for mol_str in mol_list:
                mol, r = mol_str.split(":")
                if mol != molecule:
                    continue
                
                fr, to = map(int,r.split("_"))   # these are sequence indices
                seq3_idx = 0     
                for idx in range(fr, to + 1):                    
                    row_idx = math.floor((idx-1)/nc)
                    col_idx = idx - 100 * row_idx - 1
                    color_matrix[row_idx, col_idx] = sequence_color(n_matches, seq3[seq3_idx])
                    seq3_idx += 1

        output_mol_file = output_file + "_" + molecule
        matrix_plot(color_matrix, row_labels, col_labels, mod_text_dict, output_mol_file)



def iTRAQ_quantitation():
    rep_masses = [pgv.nt_fragment_dict["light"][tag].mass_dict["Reporter"] for tag in pgv.iTRAQ_tag_set]

    # for ms2_key, mdict in pgv.ms2_match_dict.items():
    for ukey, udict in pgv.unpacked_match_dict.items(): # add iTRAQ to unpacked dict
        ms2_key = udict["ms2_key"]
        idict = pgv.ms2_dict[ms2_key].ms2
        key_list = list(idict.keys())
        rep_keys = [min(key_list, key=lambda x:abs(x-r)) for r in rep_masses]
#TODO put in delta for missing values
        rep_intensities = [idict[key] for  key in rep_keys]
        sum_rep_I = sum(rep_intensities)
        frac_rep_I = [round(x/sum_rep_I,3) for x in rep_intensities]
        for tag, fraction in zip(pgv.iTRAQ_tag_set, frac_rep_I):
            udict[tag] = fraction
        udict["sum_iTRAQ"] = sum_rep_I

def plot_ms2_spectra():
    timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M")
    plot_folder = os.path.join(pgv.job_dir,"spectrum_plots_" + timestamp)
    
    if not os.path.isdir(plot_folder):
        os.makedirs(plot_folder)
    
    ctr = 0
    for mkey, mdict in pgv.match_dict.items():
        if ctr > 0:
            break
        ctr += 1
        ukey = mdict["ukey"]
        print("ms_plot: ", ukey, generate_mod_seq(pgv.unpacked_match_dict[ukey]["frag3"]))
        # ms2_plot(plot_folder, pgv.unpacked_match_dict[ukey])
        ms2_plot(plot_folder, ukey)
    
    try: # cleanup
        os.remove(os.path.join(pgv.job_dir, 'draft.png'))
    except OSError:
        pass

####################

def ppm_offset_plot():
    
    ms1_ppm_list = np.asarray([ppm_offset(mdict["mz1"], mdict["mz_exp"]) for key, mdict in pgv.top_match_dict.items()])
    ms2_ppm_list = np.asarray([ppm_offset(idict["theo_mz2"], idict["obs_mz2"])                                 
                    for key, mdict in pgv.top_match_dict.items()
                    for idx, idict in mdict["matched_ions"].items()])
    
    ms1_mean = round(ms1_ppm_list.mean(),2)
    ms1_sd = round(ms1_ppm_list.std(),2)
    ms2_mean = round(ms2_ppm_list.mean(),2)
    ms2_sd = round(ms2_ppm_list.std(),2)
    
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 10
    xsize = 8  # spectrum plot size in inches
    ysize = 11    
    fig, ax = plt.subplots(2,1, figsize = (xsize, ysize))
    
    ax[0].hist(ms1_ppm_list,bins=30, range = (-pgv.MS1_ppm, pgv.MS1_ppm), edgecolor='black', linewidth=1.2)
    ax[1].hist(ms2_ppm_list,bins=30, range = (-pgv.MS2_ppm, pgv.MS2_ppm), edgecolor='black', linewidth=1.2)


    fig.suptitle("MS1 and MS2 match ppm offset distributions " + pgv.job_dir.split("/")[-1])
    plot_filename = "offset_distribution_plot_" + pgv.job_dir.split("_")[-1] + ".pdf"
    plot_folder = os.path.join(pgv.job_dir,plot_filename)

    
    ax[0].set_title("MS1 offsets: " + str(ms1_mean) + " +/- " + str(ms1_sd) + " ppm")
    ax[1].set_title("MS2 offsets: " + str(ms2_mean) + " +/- " + str(ms2_sd) + " ppm")
    ax[0].set_xlim(-pgv.MS1_ppm, pgv.MS1_ppm)
    ax[1].set_xlim(-pgv.MS2_ppm, pgv.MS2_ppm)
    
    plt.savefig(plot_folder)   
    plt.close(fig)
 

def find_top_target_decoy():
    top_targets_decoys = {}
    for ms2_key, ms2_dict in pgv.ms2_match_dict.items():
        if len(ms2_dict.keys()) == 0:
            continue
        target = {"rank": 999, "Sp": 0.0}
        decoy = {"rank": 999, "Sp": 0.0}
        for idx, mdict in ms2_dict.items():
            Sp = mdict["Sp"]
            rank = mdict["Sp_rank"]
            ft = mdict["frag_type"]
            if ft == "target":
                if rank < target["rank"]:
                    target = {"rank": rank, "Sp": Sp}
            elif ft == "decoy":
                if rank < decoy["rank"]:
                    decoy = {"rank": rank, "Sp": Sp}
            else:
                pass
        top_targets_decoys[ms2_key] = {"target": target, "decoy": decoy}

    return top_targets_decoys        
            
def top_Sp_histogram():
    
    Sp_dict = find_top_target_decoy()
    
    target_Sp_list = np.asarray([sdict["target"]["Sp"] for sdict in Sp_dict.values()])
    decoy_Sp_list = np.asarray([sdict["decoy"]["Sp"] for sdict in Sp_dict.values()])
        
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 10
    xsize = 8  # spectrum plot size in inches
    ysize = 11    
    fig, ax = plt.subplots(2,1, figsize = (xsize, ysize))
    
    max_Sp = max([max(target_Sp_list),max(decoy_Sp_list)])
    ax[0].hist(target_Sp_list,bins=30, range = (0, max_Sp), color = "blue", 
               alpha = 0.5, edgecolor='black', linewidth=1.2, label = "Targets")
    ax[0].hist(decoy_Sp_list,bins=30, range = (0, max_Sp), color = "red", 
               alpha = 0.5, edgecolor='black', linewidth=1.2, label = "Decoys")
    ax[0].legend()
    fig.suptitle("Top Target and Decoy Sp distributions " + pgv.job_dir.split("/")[-1])
    plot_filename = "Top_Sp_distribution_plot_" + pgv.job_dir.split("_")[-1] + ".pdf"
    plot_folder = os.path.join(pgv.job_dir,plot_filename)
    
    for idx, sdict in Sp_dict.items():
        if sdict["target"]["Sp"] > sdict["decoy"]["Sp"]:
            ax[1].plot([sdict["target"]["Sp"]],sdict["decoy"]["Sp"] , "ob")
        else:
            ax[1].plot([sdict["target"]["Sp"]],sdict["decoy"]["Sp"] , "or")

    blue_point = Line2D([0], [0], label='Target Top Sp', marker='o', 
              markerfacecolor='b', linestyle='')
    red_point = Line2D([0], [0], label='Decoy Top Sp', marker='o', 
              markerfacecolor='r',linestyle='')
       
    plt.legend(handles=[blue_point, red_point])

    ax[0].set_xlim(0, max_Sp)
    ax[1].set_xlim(0, max_Sp)
    ax[1].set_ylim(0, max_Sp)
    ax[1].set_xlabel("Target Sp")
    ax[1].set_ylabel("Decoy Sp")
    
    plt.savefig(plot_folder)   
    plt.close(fig)                
