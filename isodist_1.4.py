#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 14:31:53 2022

@author: jrwill
"""

# python version of isodist_3.2b   march 2017 version
import sys
import os
from scipy.fft import rfft, irfft,rfftfreq
from scipy.optimize import least_squares
import numpy as np
from matplotlib import pyplot as plt
import copy
import time

package_dir = "/Users/jrwill/prog/Pytheas_Folder/Pytheas_Qt_root/pytheas_dev/"
pytheas_dir =  "/Users/jrwill/prog/Pytheas_Folder/Pytheas_Qt_root/Pytheas_Qt/"
sys.path.insert(0, package_dir)
sys.path.insert(0, pytheas_dir)

from isodist_global_variables import igv
from isodist_functions import (create_atom_dict, create_model,
                               create_atom_mu, initialize_residue_mu, calc_isodist, 
                               calculate_residuals, calculate_mol_form, 
                               plot_obs_calc)
from minispectrum_functions import build_minispectra, plot_minispectrum
from mod_seq_functions import parse_mod_seq
from pytheas_global_vars import pgv,pgc
from pytheas_objects import fragment_sequence

# read Mzml file and extract minispectra for fitting

pythid_file = "/Users/jrwill/prog/Pytheas_Folder/Pytheas_data/human_28S_RNA_for_massacre/Match_Job_099/HEK293T_28S_meth100_T1_P1-A-2_01_14293_massacre.csv"
ms_data_file = "/Users/jrwill/prog/Massacre_backup/Massacre_RNA_test_set/HEK293T_28S_meth100_T1_P1-A-2_01_14293.mzML"

match_df, minispec = build_minispectra(pythid_file, ms_data_file)


#TODO QC on RT peak fits   45 /1440 don't fit
for i in range(len(minispec)):
    fit_rt_peak(i)
    
for i in range(len(minispec)):
    sum_spectra(i)

    
    
    

#TODO make this a function
# determine active residues

#TODO fix m66A

frag3_list = []
active_residues = []
for idx, row in match_df.iterrows():
    end5 = row["5'-end"]
    seq =  row["Pytheas_sequence"].replace("m66A", "m6,6A")
    end3 =  row["3'-end"]
    seq3 = parse_mod_seq(seq)
    frag3 = [pgv.end_dict["end5"][end5]] + seq3 + [pgv.end_dict["end3"][end3]]
    frag3_list.append(frag3)
    for res in frag3:
        if res not in active_residues:
            active_residues.append(res)

print("active residues: ", active_residues)
            
            





# infile = "/Users/jrwill/prog/Massacre_backup/massacre_src/bleedingEdge/sample_data/15N_pulse.in"
opts = ["showfit", "showguess"]
# cwd = '/Users/jrwill/prog/Massacre_backup/massacre_src/bleedingEdge/sample_data'



# igv.batchfile = "15N_pulse.batch" # batch input file of peaks to be fit
igv.atomfile = "pytheas_atom_definitions.txt" # atom definition file  = atom_definitions.txt
igv.model_file = "pytheas_label_model.txt"
# igv.resfile = "res_frac_15N_def.txt"   # residue definition file

igv.batchfile = "test_15N_pulse.txt"
print("batchfile ", igv.batchfile)    
igv.batchout = igv.batchfile.split(".")[0] + ".csv"    # for output of results
igv.csvfile = open(igv.batchout,"w")

# mz = np.zeros(0)
# for i in range(igc.npt):            # mz axis
#     mz = np.append(mz,float(i+1)/igc.scale_mz)

#TODO check if 1st point should be zero or, as above

igv.mz = np.linspace(0, igv.npt/igv.scale_mz, igv.npt)  # calc mz axis
igv.cmz = rfftfreq(igv.npt) # complex calc mz axis

create_atom_dict()  # input atom parameters   jgv.atom_dict
create_model() # igv.species_dict  igv.atom_order_dict igv.atom_fix_dict igv.model

# create_residue_dicts()  # set up residue library and build model  DONT NEED THIS...using nt_mod_defs 
create_atom_mu()  # generate mu domain spectra for atoms  (natom_types, ncp) 
initialize_residue_mu(active_residues)  # generate mu domain spectra for residues



#TODO make function to assemble fit from model
#     return parlabel, x_bounds, nfitpar

def initialize_fit_model():
    igv.parlabel = ["B","OFF","GW"]
    igv.x_init = [1.0, 0.01, 0.003]  # initial values for B, OFF, GW
    lower_bounds = [-1000.0,-1.0,0.001]
    upper_bounds = [np.inf,1.0,0.03]
    # nfitpar = 2   # why is this not 3???  perhaps this is index in python
    fitpar_idx = 2
    
    # amplitude parameters
    alab = "AMP_"
    igv.model.amp_idx = []
    for i in range(igv.model.n_species):
        # nfitpar = nfitpar + 1
        fitpar_idx += 1
        igv.model.amp_idx.append(fitpar_idx)
        igv.parlabel.append(alab + igv.model.species_lab[i])
        igv.x_init.append(igv.model.amp[i])
        lower_bounds.append(0.0)
        upper_bounds.append(np.inf)
    
    # fractional atom parameters
    flab = "FRC_"
    igv.model.frac_idx = []
    for i in range(igv.model.n_var_atom):
        # nfitpar = nfitpar + 1
        fitpar_idx += 1
        igv.model.frac_idx.append(fitpar_idx)
        igv.parlabel.append(flab + igv.model.var_atom[i])
        fiso_list = igv.atom_dict[igv.model.var_atom[i]]["fiso"]
        igv.x_init.append(igv.model.frac[i])
        lower_bounds.append(0.0)
        upper_bounds.append(1.0)
    
    igv.xbounds = (lower_bounds,upper_bounds)
 
    print("number of params to fit :",len(igv.x_init),igv.x_init)
    
    
    # build header for csv output
    igv.csvlabels = igv.outlabels + igv.parlabel
    
    igv.csvlabels = igv.csvlabels + ["max_fit","min_fit","max_rsd","min_rsd"]
    nfit = 11  
    nwr = 23
    for nf in range(nfit):
        igv.csvlabels.append("avg_fit" + str(nf+1))
    for nw in range(nwr):
        igv.csvlabels.append("avg_wr" + str(nw+1))
        
    csvstr = ",".join(igv.csvlabels)
    igv.csvfile.write(csvstr + '\n')  # open file for write


def fit_minispectrum(match_df_idx):
    match_row = match_df.iloc[match_df_idx]
    ms = minispec[match_df_idx]
    
    end5 = match_row["5'-end"]
    seq =  match_row["Pytheas_sequence"]
    end3 =  match_row["3'-end"]
    seq3 = parse_mod_seq(seq)
    z = abs(ms.charge)
    moz = match_row["m/z"]

    frag3 = [pgv.end_dict["end5"][end5]] + seq3 + [pgv.end_dict["end3"][end3]]
        
    print("SEQUENCE, z = ", frag3, z)
    
    peakfile = "_".join([str(match_row["pythid"]), ms.molecule_ID, ms.seq, str(z)]) + ".txt"
    seqlab = ms.molecule_ID

    mz_hd = z * ms.mzlo
    print("mz_hd = ",mz_hd)
    yobs = ms.rt_slice
    xobs = ms.mini_mz
    nmz = len(yobs)
    
    # mz_ptr has indices of experimental points:  calc spectrum is hi res
    
    ri_array = ms.rt_slice  # experimental spectrum
    rmz_array = ms.mini_mz

    nmz = 0
    mz_ptr = np.zeros(len(rmz_array), dtype = int)
    xobs = np.zeros(len(rmz_array))
    yobs = np.zeros(len(rmz_array))
    mz_hd = rmz_array[0] * z
    
    for rmz, ri in zip(rmz_array, ri_array):
        n = int((rmz * z - mz_hd)*igv.scale_mz + 0.5)

        mz_ptr[nmz] = n
        xobs[nmz] = (rmz * z - mz_hd) * igv.scale_mz
        yobs[nmz] = ri
        nmz += 1
        
    igv.mz_ptr = mz_ptr
    
    x_init = igv.x_init.copy()
    x_init[3] = max(yobs)
    x_init[4] = max(yobs)
    
#TODO are GW and AMP anticorrelated?  
    
    if "showguess" in opts:     # calc initial guess
        
        mz_spec, spec = calc_isodist(x_init, "spectrum", frag3, yobs, mz_hd, z)
        resid = calc_isodist(igv.x_init, "residuals", frag3, 12*yobs, mz_hd, z)
        plot_obs_calc(mz_spec, yobs, mz_spec, spec, seq + " Initial guess")
     
    prelim_end = time.time()
    fit_start = prelim_end

    lsq_soln = least_squares(calc_isodist, igv.x_init, verbose = 2, bounds = igv.xbounds, 
                             x_scale = 'jac', max_nfev=100,
                             args =("residuals", frag3, yobs, mz_hd, z))

    x_fit = lsq_soln.x
    chisq = 2.0 *lsq_soln.cost/(igv.sig_global*igv.sig_global*nmz)  # factor of 2 to make cost correpsond to chisquared
        
    xfit = [round(x,6) for x in x_fit]
    
    fit_end = time.time()
    print("fit time = ", fit_end - fit_start)
    print(" solution = ", xfit)
           
    mz_fit, fit = calc_isodist(x_fit, "spectrum", frag3, yobs, mz_hd, z)
    resid = calc_isodist(igv.x_init, "residuals", frag3, yobs, mz_hd, z)
    # print("resid type = ",type(resid))
    if "showfit" in opts:
        plot_obs_calc(mz_fit, yobs, mz_fit, fit, seq + " Final Fit")
   
    max_fit, min_fit, max_rsd, min_rsd, avg_fit, avg_wr = calculate_residuals(x_fit, fit, resid)
        
    # add to csv batchout
    csvdata = [peakfile, seqlab, seq, moz, z, chisq, mz_hd]
    print(type(csvdata),type(x_fit))
    csvdata = csvdata + x_fit.tolist() + [max_fit, min_fit, max_rsd, min_rsd] + avg_fit + avg_wr
    
    csvstr = []
    for csv in csvdata:
        csvstr.append(str(csv))
  
    igv.csvfile.write(",".join(csvstr) + '\n')

# todo store this in minispectrum
# write out fit file
    fitfile = open(peakfile.replace(".txt",".fit"),"w")
    for i in range(len(mz_fit)):
        fitfile.write(str(mz_fit[i])+","+str(fit[i]) + '\n')
    fitfile.close()
    
# write out fit plot
    plotfile = peakfile.replace(".txt",".png")
    plt.plot(mz_fit, yobs,".k")
    plt.ylabel("amplitude")
    plt.xlabel("m/z")
    plt.title(seq + " m0 = " + str(round(mz_hd,3)))
    plt.plot(mz_fit, fit, "-r")
    plt.savefig(plotfile)
    plt.clf()

    return mz_fit, fit

initialize_fit_model()

for i in range(len(match_df)):
    print()
    print("############## fitting # ", i)
    fit_minispectrum(i)

igv.csvfile.close()

