#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 13:46:53 2024

@author: jrwill
"""
import numpy as np

from pytheas_global_vars import pgv, pgc


def ppm_range(value, ppm): # Calculate incertitude on MS1/MS2 masses equivalent to given ppm
    return ppm / 1000000 * value

def ppm_offset(measured_mass, theoretical_mass): # Calculate the ppm offset between matching m/z values for MS1/MS2 ions
    difference = np.float64(measured_mass) - np.float64(theoretical_mass)
    ppm_offset = difference / np.float64(theoretical_mass) * 1000000

    return round(ppm_offset, 1)

def find_losses_freebases(MS2_ions):  # NOT USED
    """
    Finds the M-xx ions and free bases among the MS2 matches ions per each precursor ion
    """
    output_loss_masses = []
    for x in MS2_ions:
        if 'M-' in x[1].split("(")[0] or len(x[1].split("(")[0]) == 1:
            output_loss_masses.append(np.float64(x[2]))

    return output_loss_masses


def consecutive_series(idict, ion_series):
    """
    Return beta increment for a particular ion series
    beta is determined as:
            if no consecutive matches: beta = 0
            if at least one consecutive match: beta = sum(beta_increment + (alpha * consecutive_matches ) *
            beta_increment) over all consecutive matches in the same ion series
    """
    beta, consec_flag = 0, 0
    
    # list_MS2 = ["".join(key.split(":")[0:2]) for key in idict.keys()]  #series:index

    list_MS2 = [idict[key]["index"] for key in idict.keys() if idict[key]["series"] == ion_series]
    
    for i in range(1, 20): # should change this to max sequence length in pgv
    
         # Add the support for a-b, y-P and z-P that have the number after the first character and not at the end
         # if len(ion_series) == 1:
         #     current, consecutive = ion_series + str(i), ion_series + str(i + 1)
         #     consec_flag += 1
         # else:  # for y-P etc...
         #     current, consecutive = ion_series[0] + str(i) + ion_series[1:], ion_series[0] + str(i + 1) + ion_series[1:]
         current, consecutive = i, i+1 
         if current in list_MS2 and consecutive in list_MS2:
             beta += pgv.beta_increment + (pgv.alpha * consec_flag) * pgv.beta_increment
             consec_flag += 1
    
         else:
             consec_flag = 0
    # print("consec_series_new: ", ion_series, beta)     
    if pgv.weight_beta == 'y':
        beta = beta * pgc.iw_dict[ion_series]        
    return beta


def sumI(ion_dict): # ion_dict is indexed dict of <ms2_match> objects
    sum_int = 0.0
    mgf_peaks = []
    
    # for ms2_ion, idict in ion_dict.items():
    #     print("sumI:", ms2_ion, idict, type(idict))
    
    # print("sumI ion_dict keys", ion_dict.keys())
    for ms2_ion, idict in ion_dict.items():
        # print("ms2_ion", ms2_ion)
        # print("SumI", ms2_ion, idict.__dict__)
        # ms2_match = np.float64(idict["ms2_mass"])
        if "obs_mz2" in idict:
            ms2_match = idict["obs_mz2"]
        else:
            ms2_match = idict["mz2"]

        if float(pgv.MS2_mzlow) < ms2_match < float(pgv.MS2_mzhigh): # this should be done on spectrum object
            # series = ms2_ion.split("_")[0]
            if ms2_match not in mgf_peaks: # sum only once
                mgf_peaks.append(ms2_match)

# TODO free base and M- exclusion from sumI  
                wgt = 1.0
                if pgv.weight_ions == "y":
                    if idict["series"] in pgc.iw_dict:
                        wgt = pgc.iw_dict[idict["series"]]
                
                if "obs_int" in idict:
                    sum_int += idict["obs_int"] * wgt
                else:
                    sum_int += wgt

    return sum_int
                
def sumI_all(ms2_spectrum, thresh):
    sum_int = 0.0
    mgf_peaks = []
    ion_dict = ms2_spectrum.ms2
    max_int = float(ms2_spectrum.max_int)
    # print("SUMIALL: len ion dict", len(ion_dict.keys()))
    
    imax = 0.0
    for mz, intensity in ion_dict.items():
        if intensity > imax:
            imax = intensity
    
    for mz, intensity in ion_dict.items():
        if float(pgv.MS2_mzlow) < mz < float(pgv.MS2_mzhigh) and intensity > thresh:
            # this should be done in spectrum object
            if mz not in mgf_peaks: # sum only once
                mgf_peaks.append(mz)
                if ion_dict[mz] > imax:
                    imax = ion_dict[mz]
    
    if len(mgf_peaks) == 0:
        imax = 1.0

    sum_int = sum([ion_dict[x] for x in mgf_peaks])

    return sum_int, max_int

def n_calc(ion_dict, mod_bases, thresh):
    """
    Return the value n (number of MS2 ions) for the scoring function disregarding multiple charged species for the same series
    (e.g. a2(-1) and a2(-2) will count as 1)
    """
    n, unique_series = 0, []

    # MS2 ions following specific criteria are recursively added to the n output
    for ms2_ion, idict in ion_dict.items():

        # if pgv.MS2_mzlow < np.float64(l[0]) < pgv.MS2_mzhigh:

        # series = l[1].split('(')[0]
        # ser,idx,z = ms2_ion.split(":")
        # series = ser + idx
        series = idict["series"] + str(idict["index"])

        # Exclude all ions from neutral/charged losses and free bases
        if 'M-' not in series and series != 'G' and series != 'A' and series != 'C' and series != 'U':

            # Exclude free bases from modified nucleotides
            if mod_bases:
                if series not in mod_bases and series not in unique_series:
                    n += 1
                    unique_series.append(series)
            else:
                if series not in unique_series:
                    n += 1
                    unique_series.append(series)

    return n

def L_calc(ion_dict, mod_bases, thresh):
    """
    Return the value n (number of MS2 ions) for the scoring function disregarding multiple charged species for the same series
    (e.g. a2(-1) and a2(-2) will count as 1)
    """
    n, unique_series = 0, []

    imax = 0.0
    # for ms2_ion, idict in ion_dict.items():
        # if idict.intensity > imax:
        #     imax = idict.intensity

    # for ms2_ion, idict in ion_dict.items():
    #     if idict["obs_int"] > imax:
    #         imax = idict["obs_int"]

    # MS2 ions following specific criteria are recursively added to the n output
    for ms2_ion, idict in ion_dict.items():

        # if pgv.MS2_mzlow < np.float64(l[0]) < pgv.MS2_mzhigh:

        # series = l[1].split('(')[0]
        # ser,idx,z = ms2_ion.split(":")
        # series = ser + idx
        series = idict["series"] + str(idict["index"])

        # Exclude all ions from neutral/charged losses and free bases
        # print("L_calc", idict, imax, thresh)
        # if 100.0 * idict["obs_int"]/imax < thresh:
        #     continue
        
        if 'M-' not in series and series != 'G' and series != 'A' and series != 'C' and series != 'U':

            # Exclude free bases from modified nucleotides
            if mod_bases:
                if series not in mod_bases and series not in unique_series:
                    n += 1
                    unique_series.append(series)
            else:
                if series not in unique_series:
                    n += 1
                    unique_series.append(series)

    return n
